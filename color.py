# This file provides utility functions for constructing matplotlib color maps that are adapted for use in plotting modular forms

from matplotlib.colors import LinearSegmentedColormap, to_rgb, to_hex

def smooth_gradient(*colors, name=None):
    """
    A color map using the given colors without sharp boundaries.

    INPUT:

    - ``colors`` -- colors provided either via hex values or rgb tuples

    OUTPUT:

    A matplotlib colormap object
    """
    colors = [to_rgb(color) for color in colors]
    if len(set(colors)) < 2:
        raise ValueError("You must provide at least two colors")
    if name is None:
        name = "".join(to_hex).replace("#","")
    n = len(colors)
    colors = colors + [colors[0]]
    breaks = [float(i/n) for i in range(n+1)]
    return LinearSegmentedColormap(name, segmentdata={
        "red": [(x, col[0], col[0]) for (x,col) in zip(breaks, colors)],
        "green": [(x, col[1], col[1]) for (x,col) in zip(breaks, colors)],
        "blue": [(x, col[2], col[2]) for (x, col) in zip(breaks, colors)]})

gradients = {
    "MIT": smooth_gradient("#750014", "#8b959e", name="MITGrad"),
}
