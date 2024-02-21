
import os
from color import gradients
from matplotlib import colormaps
from modform import ModularForm
from hyperbolic import hyperbolic_model
from sage.plot.complex_plot import complex_plot
from sage.rings.infinity import infinity

def plot_modform(label, resolution=500, cmap=None, hmodel=None):
    if hmodel is None:
        hmodel = hyperbolic_model("poincare")
    else:
        hmodel = hyperbolic_model(hmodel)
    if cmap is None:
        cmap = gradients["MIT"]
    elif cmap in gradients:
        cmap = gradients[cmap]
    elif cmap in colormaps:
        cmap = colormaps[cmap]
    modform = ModularForm(label)
    def evaluate(z):
        if (hmodel.disk and abs(z) >= .995) or (not hmodel.disk and z.imag() < 0.005):
            return infinity
        return modform.evaluate(hmodel.convert(z))

    P = complex_plot(
        evaluate,
        hmodel.xrange,
        hmodel.yrange,
        cmap=cmap,
        dark_rate=0.00,
        contoured=True,
        plot_points=resolution,
        aspect_ratio=1)
    P.axes(False)
    return P

def save_modform_plot(label, resolution=500, cmap=None, hmodel=None, filename=None):
    if filename is None:
        filename = os.path.expanduser(f"~/Downloads/{label}_{resolution}.png")
    P = plot_modform(label, resolution=resolution, cmap=cmap, hmodel=hmodel)
    P.save(filename, figsize=[20,20], axes=False)
