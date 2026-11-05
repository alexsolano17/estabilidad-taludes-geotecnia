"""Línea de comandos.

    python -m estabilidad_taludes analizar casos/talud_cobertura.json
    python -m estabilidad_taludes diseno casos/talud_cobertura.json --desde 30 --hasta 55 --paso 2.5
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from . import __version__
from .analisis import ResultadoAnalisis, analizar, angulo_recomendado, barrido_angulos
from .casos import cargar_caso
from .graficos import grafico_convergencia, grafico_diseno, grafico_histograma, grafico_talud


def _redondear(obj, cifras: int = 6):
    if isinstance(obj, float):
        return float(f"{obj:.{cifras}g}")
    if isinstance(obj, dict):
        return {k: _redondear(v, cifras) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_redondear(v, cifras) for v in obj]
    return obj


def resumen(r: ResultadoAnalisis) -> dict:
    """Resultados principales en un diccionario serializable (sin tiempos, para que sea reproducible)."""
    caso = r.caso
    datos = {
        "version": __version__,
        "caso": caso.nombre,
        "descripcion": caso.descripcion,
        "geometria": {**caso.geometria, **({"angulo": r.angulo} if r.angulo is not None else {})},
        "material": asdict(caso.material),
        "nivel_freatico": caso.nivel_freatico,
        "determinista": {
            "fs_bishop": r.bishop.critica.fs,
            "fs_fellenius": r.fellenius.critica.fs,
            "fs_bishop_sin_agua": r.fs_seco,
            "superficie_critica_bishop": asdict(r.bishop.critica),
            "superficie_critica_fellenius": asdict(r.fellenius.critica),
            "dovelas": r.bishop.n_dovelas,
            "circulos_malla": r.bishop.evaluadas,
            "circulos_validos": r.bishop.validas,
        },
    }
    if r.monte_carlo is not None:
        datos["probabilistico"] = {
            "simulaciones": int(r.monte_carlo.c.size),
            "semilla": caso.semilla,
            "variables": [asdict(v) for v in caso.variables],
            "correlacion_c_phi": caso.rho_c_phi,
            "superficies_evaluadas": len(r.monte_carlo.superficies),
            "avisos_bishop": r.monte_carlo.avisos_bishop,
            "valores_recortados": r.monte_carlo.recortadas,
            "minimo_entre_superficies": asdict(r.est_min),
            "superficie_critica_fija": asdict(r.est_fija),
            "rosenblueth_superficie_fija": {
                "media": r.rosenblueth.media,
                "desviacion": r.rosenblueth.desviacion,
                "puntos": [
                    {**v, "peso": float(p), "fs": float(f)}
                    for v, p, f in zip(r.rosenblueth.valores, r.rosenblueth.pesos, r.rosenblueth.fs)
                ],
            },
        }
    if caso.criterios is not None:
        datos["criterios"] = asdict(caso.criterios)
        datos["cumple"] = r.cumple()
    return _redondear(datos)


def _pct(p: float) -> str:
    return f"{100 * p:.2f} %"


def resumen_markdown(r: ResultadoAnalisis) -> str:
    c, m = r.caso, r.caso.material
    filas = [
        f"# Resultados: {c.nombre}",
        "",
        c.descripcion,
        "",
        "## Datos",
        "",
        "| Parámetro | Valor |",
        "|---|---|",
    ]
    g = c.geometria
    if c.es_simple:
        filas += [f"| Altura | {g['altura']:g} m |", f"| Ángulo de la cara | {r.angulo:g}° |"]
    filas += [
        f"| c' | {m.c:g} kPa |",
        f"| φ' | {m.phi:g}° |",
        f"| γ / γsat | {m.gamma:g} / {m.gamma_saturado:g} kN/m³ |",
        f"| Nivel freático | {'sí' if r.talud.nivel_freatico is not None else 'no'} |",
    ]
    for v in c.variables:
        nombre, media = ("c'", f"{v.media:g} kPa") if v.nombre == "c" else ("φ'", f"{v.media:g}°")
        filas.append(f"| {nombre} aleatoria | {v.distribucion}, media {media}, COV {v.cov:.0%} |")
    if c.variables:
        filas.append(f"| Correlación ρ(c', φ') | {c.rho_c_phi:g} |")
    s = r.bishop.critica
    filas += [
        "",
        "## Análisis determinista",
        "",
        "| Resultado | Valor |",
        "|---|---|",
        f"| FS Bishop simplificado | **{s.fs:.3f}** |",
        f"| FS Fellenius (ordinario) | {r.fellenius.critica.fs:.3f} |",
    ]
    if r.fs_seco is not None:
        filas.append(f"| FS Bishop sin nivel freático | {r.fs_seco:.3f} |")
    filas += [
        f"| Círculo crítico | centro ({s.xc:.2f}, {s.yc:.2f}) m, radio {s.radio:.2f} m |",
        f"| Salida / entrada | x = {s.x_salida:.2f} m / x = {s.x_entrada:.2f} m |",
        f"| Círculos válidos en la malla | {r.bishop.validas:,} de {r.bishop.evaluadas:,} |".replace(",", " "),
    ]
    if r.est_min is not None:
        e, f, ros = r.est_min, r.est_fija, r.rosenblueth
        filas += [
            "",
            "## Análisis probabilístico",
            "",
            f"{e.n:,} simulaciones de Monte Carlo con semilla {c.semilla}; "
            f"{len(r.monte_carlo.superficies)} superficies evaluadas por muestra.".replace(",", " "),
            "",
            "| Resultado | Mínimo entre superficies | Superficie crítica fija |",
            "|---|---|---|",
            f"| Media del FS | {e.media:.3f} | {f.media:.3f} |",
            f"| Desviación del FS | {e.desviacion:.3f} | {f.desviacion:.3f} |",
            f"| Percentil 5 % del FS | {e.p05:.3f} | {f.p05:.3f} |",
            f"| **P(FS < 1)** | **{_pct(e.pf)}** | {_pct(f.pf)} |",
            (
                f"| IC 95 % de P(FS < 1) | {_pct(e.pf_ic95[0])} – {_pct(e.pf_ic95[1])} | "
                f"{_pct(f.pf_ic95[0])} – {_pct(f.pf_ic95[1])} |"
            ),
            f"| β (FS normal) | {e.beta_normal:.2f} | {f.beta_normal:.2f} |",
            f"| β (FS lognormal) | {e.beta_lognormal:.2f} | {f.beta_lognormal:.2f} |",
            "",
            (
                f"Control con Rosenblueth (superficie fija, {len(ros.fs)} puntos): media {ros.media:.3f}, "
                f"desviación {ros.desviacion:.3f}."
            ),
        ]
    if c.criterios is not None:
        cumple = r.cumple()
        ok = {True: "cumple", False: "no cumple", None: "—"}
        filas += [
            "",
            "## Criterios de aceptación",
            "",
            c.criterios.descripcion,
            "",
            "| Criterio | Límite | Resultado | Estado |",
            "|---|---|---|---|",
            f"| FS mínimo | ≥ {c.criterios.fs_min:g} | {s.fs:.3f} | {ok[cumple['fs']]} |",
        ]
        if r.pf is not None:
            filas.append(f"| P(FS < 1) máxima | ≤ {_pct(c.criterios.pf_max)} | {_pct(r.pf)} | {ok[cumple['pf']]} |")
    return "\n".join(filas) + "\n"


def guardar_resultados(r: ResultadoAnalisis, carpeta: Path) -> list[Path]:
    carpeta.mkdir(parents=True, exist_ok=True)
    rutas = [grafico_talud(r, carpeta / "talud_critico.png")]
    (carpeta / "resumen.json").write_text(json.dumps(resumen(r), ensure_ascii=False, indent=2), encoding="utf-8")
    (carpeta / "resumen.md").write_text(resumen_markdown(r), encoding="utf-8")
    rutas += [carpeta / "resumen.json", carpeta / "resumen.md"]
    if r.monte_carlo is not None:
        rutas.append(grafico_histograma(r, carpeta / "histograma_fs.png"))
        rutas.append(grafico_convergencia(r, carpeta / "convergencia_pf.png"))
        mc = r.monte_carlo
        with open(carpeta / "muestras.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["c_kPa", "phi_grados", "fs_superficie_critica", "fs_minimo", "superficie_del_minimo"])
            for fila in zip(mc.c, mc.phi, mc.fs_fija, mc.fs_min, mc.superficie_min):
                w.writerow([f"{fila[0]:.3f}", f"{fila[1]:.3f}", f"{fila[2]:.4f}", f"{fila[3]:.4f}", int(fila[4])])
        rutas.append(carpeta / "muestras.csv")
    return rutas


def _imprimir(r: ResultadoAnalisis) -> None:
    print(f"Caso: {r.caso.nombre}")
    print(f"  FS Bishop = {r.fs:.3f}   FS Fellenius = {r.fellenius.critica.fs:.3f}", end="")
    print("" if r.fs_seco is None else f"   FS sin agua = {r.fs_seco:.3f}")
    if r.est_min is not None:
        e = r.est_min
        print(f"  P(FS < 1) = {_pct(e.pf)} (IC 95 %: {_pct(e.pf_ic95[0])} – {_pct(e.pf_ic95[1])}), "
              f"β lognormal = {e.beta_lognormal:.2f}")
    cumple = r.cumple()
    if cumple is not None:
        print(f"  Criterios: {'CUMPLE' if cumple['global'] else 'NO CUMPLE'} ({r.caso.criterios.descripcion})")
    print(f"  Tiempo: {r.segundos:.1f} s")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m estabilidad_taludes",
        description="Factor de seguridad (Bishop simplificado) y probabilidad de falla (Monte Carlo) de taludes.",
    )
    sub = parser.add_subparsers(dest="comando", required=True)
    pa = sub.add_parser("analizar", help="analiza un caso y guarda figuras y resúmenes")
    pa.add_argument("caso", type=Path, help="archivo JSON del caso")
    pa.add_argument("--salida", type=Path, help="carpeta de resultados (por defecto resultados/<caso>)")
    pa.add_argument("--simulaciones", type=int, help="reemplaza el número de simulaciones del caso")
    pd = sub.add_parser("diseno", help="barrido del ángulo de la cara para elegir el ángulo de diseño")
    pd.add_argument("caso", type=Path)
    pd.add_argument("--desde", type=float, default=30.0)
    pd.add_argument("--hasta", type=float, default=55.0)
    pd.add_argument("--paso", type=float, default=2.5)
    pd.add_argument("--salida", type=Path)
    pd.add_argument("--simulaciones", type=int)
    args = parser.parse_args(argv)

    caso = cargar_caso(args.caso)
    carpeta = args.salida or Path("resultados") / args.caso.stem

    if args.comando == "analizar":
        r = analizar(caso, simulaciones=args.simulaciones)
        _imprimir(r)
        for ruta in guardar_resultados(r, carpeta):
            print(f"  -> {ruta}")
        return 0

    angulos = np.round(np.arange(args.desde, args.hasta + 1e-9, args.paso), 6)
    filas = barrido_angulos(caso, angulos, args.simulaciones)
    recomendado = angulo_recomendado(filas)
    carpeta.mkdir(parents=True, exist_ok=True)
    with open(carpeta / "diseno.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["angulo", "fs_bishop", "fs_bishop_sin_agua", "pf", "pf_ic95_inf", "pf_ic95_sup", "cumple_fs", "cumple_pf"])
        for fila in filas:
            ic = fila.pf_ic95 or (None, None)
            w.writerow([
                f"{fila.angulo:g}", f"{fila.fs:.4f}",
                "" if fila.fs_seco is None else f"{fila.fs_seco:.4f}",
                "" if fila.pf is None else f"{fila.pf:.5f}",
                "" if ic[0] is None else f"{ic[0]:.5f}", "" if ic[1] is None else f"{ic[1]:.5f}",
                fila.cumple_fs, fila.cumple_pf,
            ])
    titulo = f"{caso.nombre}: ángulo de diseño"
    grafico_diseno(filas, caso.criterios, recomendado, carpeta / "diseno_angulo.png", titulo)
    print(f"Caso: {caso.nombre}")
    print("  ángulo   FS     FS sin agua   P(FS<1)")
    for fila in filas:
        seco = "   —  " if fila.fs_seco is None else f"{fila.fs_seco:6.3f}"
        pf = "  —" if fila.pf is None else _pct(fila.pf)
        print(f"  {fila.angulo:5.1f}°  {fila.fs:5.3f}   {seco}       {pf}")
    print(f"  Ángulo recomendado: {'ninguno cumple' if recomendado is None else f'{recomendado:g}°'}")
    print(f"  -> {carpeta / 'diseno.csv'}\n  -> {carpeta / 'diseno_angulo.png'}")
    return 0
