"""Figuras del análisis (PNG): sección crítica, histograma del FS, convergencia y diseño."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .probabilistico import intervalo_wilson, pf_acumulada

# Paleta: superficie, tintas y serie 1-2 de la paleta categórica validada (azul, naranja).
SUPERFICIE = "#fcfcfb"
TINTA = "#0b0b0b"
TINTA_2 = "#52514e"
TENUE = "#898781"
GRILLA = "#e1e0d9"
EJE = "#c3c2b7"
AZUL = "#2a78d6"
AZUL_CLARO = "#cde2fb"
NARANJA = "#eb6834"
CRITICO = "#d03b3b"  # color de estado: sólo para la zona de falla (FS < 1)
SUELO = "#ece9e1"
DPI = 160

ESTILO = {
    "figure.facecolor": SUPERFICIE,
    "axes.facecolor": SUPERFICIE,
    "savefig.facecolor": SUPERFICIE,
    "axes.edgecolor": EJE,
    "axes.labelcolor": TINTA_2,
    "axes.titlecolor": TINTA,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRILLA,
    "grid.linewidth": 0.6,
    "grid.linestyle": "-",
    "xtick.color": TENUE,
    "ytick.color": TENUE,
    "xtick.labelcolor": TINTA_2,
    "ytick.labelcolor": TINTA_2,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "font.size": 10,
}


def _guardar(fig, ruta: Path) -> Path:
    ruta = Path(ruta)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return ruta


def _arco(sup, n: int = 200):
    x = np.linspace(sup.x_salida, sup.x_entrada, n)
    y = sup.yc - np.sqrt(np.clip(sup.radio**2 - (x - sup.xc) ** 2, 0.0, None))
    return x, y


def grafico_talud(resultado, ruta: Path, n_otras: int = 8) -> Path:
    """Sección del talud con el círculo crítico de Bishop, sus dovelas y el nivel freático."""
    talud = resultado.talud
    crit = resultado.bishop.critica
    _, _, altura = talud.cara()
    with plt.rc_context(ESTILO):
        fig, ax = plt.subplots(figsize=(9.5, 5.6))
        ax.grid(False)
        xs_otros = [crit.x_salida, crit.x_entrada]
        ys_otros = [crit.yc]
        mc = resultado.monte_carlo
        otras = [] if mc is None else mc.superficies[1 : 1 + n_otras]
        for sup in otras:
            xs_otros += [sup.x_salida, sup.x_entrada]
            ys_otros.append(sup.yc)

        x_izq = max(talud.x_min, min(xs_otros) - 0.6 * altura)
        x_der = min(talud.x_max, max(xs_otros) + 0.6 * altura)
        y_bajo = min(float(np.min(_arco(s)[1])) for s in [crit, *otras]) - 0.25 * altura
        if talud.y_estrato_firme is not None:
            y_bajo = min(y_bajo, talud.y_estrato_firme - 0.1 * altura)
        y_alto = max(np.max(talud.superficie[:, 1]) + 0.3 * altura, min(crit.yc + 0.1 * altura, 3 * altura))

        x_terreno = np.linspace(x_izq, x_der, 600)
        y_terreno = talud.y_terreno(x_terreno)
        ax.fill_between(x_terreno, y_bajo, y_terreno, color=SUELO, lw=0, zorder=1)
        ax.plot(x_terreno, y_terreno, color=TINTA_2, lw=1.6, zorder=4)
        if talud.y_estrato_firme is not None:
            ax.axhline(talud.y_estrato_firme, color=TINTA_2, lw=1.0, zorder=2)
            ax.text(x_izq, talud.y_estrato_firme, " estrato firme", ha="left", va="bottom",
                    color=TINTA_2, fontsize=8)

        for i, sup in enumerate(otras):
            xa, ya = _arco(sup)
            ax.plot(xa, ya, color=TENUE, lw=0.8, alpha=0.55, zorder=3,
                    label="Otras superficies del análisis probabilístico" if i == 0 else None)

        # Dovelas: se dibuja 1 de cada 5 para que la figura sea legible.
        bordes = np.linspace(crit.x_salida, crit.x_entrada, resultado.bishop.n_dovelas + 1)[::5]
        for xb in bordes[1:-1]:
            yb = crit.yc - np.sqrt(max(crit.radio**2 - (xb - crit.xc) ** 2, 0.0))
            ax.plot([xb, xb], [yb, talud.y_terreno(xb)], color=NARANJA, lw=0.6, alpha=0.45, zorder=5)
        xa, ya = _arco(crit)
        ax.plot(xa, ya, color=NARANJA, lw=2.4, zorder=6, label="Superficie crítica (Bishop)")
        if crit.yc <= y_alto:
            ax.plot([crit.xc], [crit.yc], marker="+", color=NARANJA, ms=10, mew=1.8, zorder=6)
            for xe in (crit.x_salida, crit.x_entrada):
                ax.plot([crit.xc, xe], [crit.yc, talud.y_terreno(xe)], color=NARANJA, lw=0.7,
                        alpha=0.6, ls=(0, (3, 3)), zorder=5)

        if talud.nivel_freatico is not None:
            y_agua = np.minimum(talud.y_freatico(x_terreno), y_terreno)
            ax.plot(x_terreno, y_agua, color=AZUL, lw=1.6, zorder=5, label="Nivel freático")
            xm = x_izq + 0.85 * (x_der - x_izq)
            ax.plot([xm], [np.minimum(talud.y_freatico(xm), talud.y_terreno(xm)) + 0.04 * altura],
                    marker="v", color=AZUL, ms=7, zorder=6)

        partes = [f"FS Bishop = {crit.fs:.3f}", f"FS Fellenius = {resultado.fellenius.critica.fs:.3f}"]
        if resultado.fs_seco is not None:
            partes.append(f"FS sin agua = {resultado.fs_seco:.3f}")
        if resultado.est_min is not None:
            partes.append(f"P(FS < 1) = {100 * resultado.est_min.pf:.2f} %")

        ax.set_xlim(x_izq, x_der)
        ax.set_ylim(y_bajo, y_alto)
        ax.set_aspect("equal")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("cota (m)")
        titulo = resultado.caso.nombre
        if resultado.angulo is not None:
            titulo += f" — cara a {resultado.angulo:g}°"
        ax.set_title(titulo, pad=24)
        ax.text(0.0, 1.015, "   ·   ".join(partes), transform=ax.transAxes, ha="left", va="bottom",
                fontsize=10, color=TINTA_2)
        ax.legend(loc="best")
        return _guardar(fig, ruta)


def grafico_histograma(resultado, ruta: Path) -> Path:
    """Distribución del FS (mínimo entre superficies candidatas) y probabilidad de falla."""
    est = resultado.est_min
    fs = resultado.monte_carlo.fs_min
    with plt.rc_context(ESTILO):
        fig, ax = plt.subplots(figsize=(8.5, 4.8))
        ax.grid(axis="x", visible=False)
        bins = np.linspace(np.percentile(fs, 0.05), np.percentile(fs, 99.95), 60)
        conteo, bordes = np.histogram(fs, bins=bins)
        densidad = conteo / (len(fs) * np.diff(bordes))
        colores = [CRITICO if b < 1.0 else AZUL for b in bordes[:-1]]
        ax.bar(bordes[:-1], densidad, width=np.diff(bordes), align="edge", color=colores,
               edgecolor=SUPERFICIE, linewidth=0.8)
        tope = 1.32 * densidad.max()
        ax.set_ylim(0, tope)
        ax.axvline(1.0, color=CRITICO, lw=1.4, ls=(0, (4, 3)))
        ax.text(1.0, 0.98 * tope, "falla: FS < 1 ", color=TINTA, ha="right", va="top", fontsize=9)
        ax.axvline(resultado.fs, color=TINTA_2, lw=1.2)
        ax.text(resultado.fs, 0.98 * tope, f" FS determinista = {resultado.fs:.3f}",
                color=TINTA_2, ha="left", va="top", fontsize=9)
        resumen = (
            f"P(FS < 1) = {100 * est.pf:.2f} %\n"
            f"IC 95 %: {100 * est.pf_ic95[0]:.2f} – {100 * est.pf_ic95[1]:.2f} %\n"
            f"media {est.media:.3f} · desviación {est.desviacion:.3f}\n"
            f"β lognormal = {est.beta_lognormal:.2f}\n"
            f"{est.n:,} simulaciones".replace(",", " ")
        )
        ax.text(0.985, 0.97, resumen, transform=ax.transAxes, ha="right", va="top", fontsize=9,
                color=TINTA, linespacing=1.5, bbox={"boxstyle": "round,pad=0.5", "fc": SUPERFICIE, "ec": GRILLA})
        ax.set_xlabel("Factor de seguridad (Bishop simplificado)")
        ax.set_ylabel("Densidad de probabilidad")
        ax.set_title("Distribución del factor de seguridad (Monte Carlo)")
        return _guardar(fig, ruta)


def grafico_convergencia(resultado, ruta: Path) -> Path:
    """Estimación de la probabilidad de falla a medida que crece el número de simulaciones."""
    fs = resultado.monte_carlo.fs_min
    n = np.arange(1, fs.size + 1)
    pf = pf_acumulada(fs)
    fallas = np.cumsum(fs < 1.0)
    puntos = np.unique(np.geomspace(10, fs.size, 300).astype(int)) - 1
    banda = np.array([intervalo_wilson(int(fallas[i]), int(n[i])) for i in puntos])
    with plt.rc_context(ESTILO):
        fig, ax = plt.subplots(figsize=(8.5, 4.4))
        ax.fill_between(n[puntos], 100 * banda[:, 0], 100 * banda[:, 1], color=AZUL_CLARO, lw=0,
                        label="Intervalo de confianza 95 % (Wilson)")
        ax.plot(n[puntos], 100 * pf[puntos], color=AZUL, lw=1.8, label="P(FS < 1) estimada")
        ax.annotate(f"{100 * pf[-1]:.2f} %", (n[-1], 100 * pf[-1]), xytext=(6, 0),
                    textcoords="offset points", va="center", color=TINTA, fontsize=9)
        ax.set_xscale("log")
        ax.set_xlim(10, fs.size * 1.6)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Número de simulaciones")
        ax.set_ylabel("Probabilidad de falla (%)")
        ax.set_title("Convergencia de la probabilidad de falla")
        ax.legend(loc="upper right")
        return _guardar(fig, ruta)


def grafico_diseno(filas, criterios, recomendado, ruta: Path, titulo: str = "") -> Path:
    """FS y probabilidad de falla según el ángulo de la cara (dos paneles con el mismo eje x)."""
    angulos = np.array([f.angulo for f in filas])
    fs = np.array([f.fs for f in filas])
    fs_seco = np.array([np.nan if f.fs_seco is None else f.fs_seco for f in filas])
    pf = np.array([np.nan if f.pf is None else 100 * f.pf for f in filas])
    with plt.rc_context(ESTILO):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 7.2), sharex=True,
                                       gridspec_kw={"height_ratios": [1.15, 1], "hspace": 0.12})
        if np.any(np.isfinite(fs_seco)):
            ax1.plot(angulos, fs_seco, color=NARANJA, lw=1.8, marker="o", ms=5, label="Sin nivel freático")
            ax1.annotate("sin agua", (angulos[0], fs_seco[0]), xytext=(-6, 0), textcoords="offset points",
                         ha="right", va="center", color=TINTA_2, fontsize=9)
        ax1.plot(angulos, fs, color=AZUL, lw=1.8, marker="o", ms=5, label="Con nivel freático")
        ax1.annotate("con agua", (angulos[0], fs[0]), xytext=(-6, 0), textcoords="offset points",
                     ha="right", va="center", color=TINTA_2, fontsize=9)
        ax2.plot(angulos, pf, color=AZUL, lw=1.8, marker="o", ms=5, label="Con nivel freático")
        if np.any(np.isfinite(pf)):
            ax2.annotate("con agua", (angulos[-1], pf[-1]), xytext=(-10, 0), textcoords="offset points",
                         ha="right", va="center", color=TINTA_2, fontsize=9)
        if criterios is not None:
            ax1.axhline(criterios.fs_min, color=TINTA_2, lw=1.1, ls=(0, (4, 3)))
            ax1.text(angulos[-1], criterios.fs_min, f"FS mínimo = {criterios.fs_min:g}", ha="right",
                     va="bottom", color=TINTA_2, fontsize=9)
            ax2.axhline(100 * criterios.pf_max, color=TINTA_2, lw=1.1, ls=(0, (4, 3)))
            ax2.text(angulos[0], 100 * criterios.pf_max, f"PF máxima = {100 * criterios.pf_max:g} %",
                     ha="left", va="bottom", color=TINTA_2, fontsize=9)
        if recomendado is not None:
            for ax in (ax1, ax2):
                ax.axvline(recomendado, color=TINTA, lw=1.0)
            ax1.text(recomendado, ax1.get_ylim()[1], f" ángulo recomendado: {recomendado:g}°",
                     ha="left", va="top", color=TINTA, fontsize=9, fontweight="bold")
        ax1.set_ylabel("FS mínimo (Bishop)")
        ax2.set_ylabel("P(FS < 1) (%)")
        ax2.set_xlabel("Ángulo de la cara del talud (°)")
        ax2.set_ylim(bottom=0)
        ax1.set_xlim(angulos[0] - 3.5, angulos[-1] + 1.0)
        ax1.legend(loc="upper right")
        ax1.set_title(titulo or "Ángulo de diseño: factor de seguridad y probabilidad de falla")
        return _guardar(fig, ruta)
