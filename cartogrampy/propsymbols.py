from __future__ import division, print_function
# Copyright. 2017. Daniel Lewis, Federica Bianco.
# https://github.com/danlewis85/cartogrampy/blob/master/LICENSE

def propsymbols(gdf, value_field, geom_field='geometry', symbol='circle', scale_factor=10, position='centroid',
                 edgecolors='k', facecolors='b', alpha=0.5):
    """Returns proportional symbols for a given GeoPandas GeoDataFrame value field as a Matplotlib PatchCollection

    Args:
        gdf (geopandas.geodataframe.GeoDataFrame): Input GeoDataFrame.
        value_field (str): Field name of values to be used to produce proportional symbols.
        geom_field (str, optional): Field name of geometry column in input GeoDataFrame, defaults to 'geometry'.
        symbol (str or int, optional): Name of proportional symbol to be created, defaults to 'circle', options available:
            'circle' - inputs accepted: 'circle' (case insensitive),'c','C', 1
            'triangle' - inputs accepted: 'triangle' (case insensitive),'t','T','Tri' (case insensitive),'3', 3
            'square' - inputs accepted: 'square' (case insensitive),'sq','s','S','4', 4
            'pentagon' - inputs accepted: 'pentagon' (case insensitive),'p','P','Pent' (case insensitive),'5', 5
            'hexagon' - inputs accepted: 'hexagon' (case insensitive),'h','H','hex' (case insensitive),'6', 6
            'heptagon' - inputs accepted: 'heptagon','septagon' (case insensitive),'hept' (case insensitive),'sept' (case insensitive),'7', 7
            'octagon' - inputs accepted: 'octagon' (case insensitive)','oct' (case insensitive),'8', 8
            'nonagon' - inputs accepted: 'nonagon' (case insensitive),'non' (case insensitive),'9', 9
            'decagon' - inputs accepted: 'decagon' (case insensitive),'dec' (case insensitive),'10',10
            'rhombus' - inputs accepted: 'rhombus' (case insensitive),'r','R'
            'lozenge' - inputs accepted: 'lozenge' (case insensitive),'l','L'
            'pentagram' - inputs accepted: 'pentagram' (case insensitive)
            'hexagram' - inputs accepted: 'hexagram' (case insensitive)
        scale_factor (numeric, optional): Arbitary Int or Float constant value to scale proportional symbols.
        position (str, optional): Center position of proportional symbols, defaults to centroid.
            'centroid' - GeoDataFrame centroid property, shape centroid.
            'center' - Centroid of GeoDataFrame envelope (bounding box) property, envelope centroid.
            'representative point' - Point returned by GeoDataFrame represntative point method, guarenteed to be inside polygon.
        edgepointlors (str, optional): Matplotlib color (defaults to black 'k')
        facecolors (str, optional): Matplotlib color (defaults to blue 'b')
        alpha (str, optional): Transparency (defaults to 0.5)

    Returns:
        Matplotlib PatchCollection of proportional Symbols
    """

    # Required matplotlib and numpy classes/functions for patch and collection
    # creation.
    from matplotlib.patches import Circle, RegularPolygon, Rectangle, Polygon
    from matplotlib.collections import PatchCollection
    from numpy import power, sin, tan, cos, pi

    # Calculate proportional symbol positions based on
    if position.lower() in ['centroid']:
        cents = gdf[geom_field].apply(
            lambda x: [x.centroid.x, x.centroid.y]).tolist()
    elif position.lower() in ['center', 'centre']:
        cents = gdf[geom_field].envelope.apply(
            lambda x: [x.centroid.x, x.centroid.y]).tolist()
    elif position.lower() in ['representative point', 'rep']:
        cents = gdf[geom_field].representative_point().apply(
            lambda x: [x.x, x.y]).tolist()
    else:
        # If position parameter not recognised, default to centroid.
        print("position parameter invalid, using 'centroid'. Options are 'centroid','center', or 'representative point'")
        cents = gdf[geom_field].apply(
            lambda x: [x.centroid.x, x.centroid.y]).tolist()

    symbol = str(symbol).lower()
    
    # Big if, elif, else statement to differentiate symbol choice.
    if symbol in ['circle', 'c', '1']:
        # Work out circular symbol radii
        radii = power(gdf[value_field] / pi, 0.5) * scale_factor
        # Create patches based on radii
        patches = [Circle(cents[i], radius=radii[i])
                   for i in range(len(radii))]

    elif symbol in ['triangle', 't', 'tri', '3']:
        # Work out (equilateral) triangle radius (circumradius)
        radii = power((gdf[value_field] * 12.0) /
                      (9.0 * power(3.0, 0.5)), 0.5) * scale_factor
        # Create patches based on radii
        patches = [RegularPolygon(cents[i], 3, radius=radii[i])
                   for i in range(len(radii))]

    elif symbol in ['square', 'sq', 's', '4']:
        # Work out square symbol side lengths
        side = power(gdf[value_field], 0.5) * scale_factor
        # Create patches based on side length
        patches = [Rectangle((cents[i][0] - side[i] / 2.0, cents[i][1] -
                              side[i] / 2.0), side[i], side[i]) 
                   for i in range(len(side))]

    elif symbol in ['pentagon', 'p', 'pent', '5']:
        # work out pentagon radius
        radii = power(
            (4.0 * tan(pi / 5.0) * gdf[value_field]) /
            5.0, 0.5) / (2.0 * sin(pi / 5.0)) * scale_factor
        # Create patches based on radii
        patches = [RegularPolygon(cents[i], 5, radius=radii[i])
                   for i in range(len(radii))]

    elif symbol in ['hexagon', 'h', 'hex','6']:
        # work out hexagon radius
        radii = power(
            (4.0 * tan(pi / 6.0) * gdf[value_field]) /
            6.0, 0.5) / (2.0 * sin(pi / 6.0)) * scale_factor
        # Create patches based on radii
        patches = [RegularPolygon(cents[i], 6, radius=radii[i])
                   for i in range(len(radii))]

    elif symbol in ['heptagon', 'septagon', 'hept', 'sept', '7']:
        # work out heptagon/septagon radius
        radii = power(
            (4.0 * tan(pi / 7.0) * gdf[value_field]) /
            7.0, 0.5) / (2.0 * sin(pi / 7.0)) * scale_factor
        # Create patches based on radii
        patches = [RegularPolygon(cents[i], 7, radius=radii[i])
                   for i in range(len(radii))]

    # NB decided not to allow 'o' as often denotes a circle
    elif symbol in ['octagon', 'oct', '8']:
        # work out octogon radius
        radii = power(
            (4.0 * tan(pi / 8.0) * gdf[value_field]) /
            8.0, 0.5) / (2.0 * sin(pi / 8.0)) * scale_factor
        # Create patches based on radii
        patches = [RegularPolygon(cents[i], 8, radius=radii[i])
                   for i in range(len(radii))]

    elif symbol in ['nonagon', 'non', '9']:
        # work out nonagon radius
        radii = power(
            (4.0 * tan(pi / 9.0) * gdf[value_field]) /
            9.0, 0.5) / (2.0 * sin(pi / 9.0)) * scale_factor
        # Create patches based on radii
        patches = [RegularPolygon(cents[i], 9, radius=radii[i])
                   for i in range(len(radii))]

    elif symbol in ['decagon', 'dec', '10']:
        # work out decagon radius
        radii = power((4.0 * tan(pi / 10.0) * gdf[value_field]) / 10.0, 0.5) / (
            2.0 * sin(pi / 10.0)) * scale_factor
        # Create patches based on radii
        patches = [RegularPolygon(cents[i], 10, radius=radii[i])
                   for i in range(len(radii))]

    elif symbol in ['rhombus', 'r']:
        # Assume diagonals are equal - is this technically a rhombus?
        radii = power(gdf[value_field] / 2.0, 0.5) * scale_factor
        # Create patches based on radii length
        patches = [RegularPolygon(cents[i], 4, radius=radii[i])
                   for i in range(len(radii))]

    elif symbol in ['lozenge', 'l']:
        # Assume acute angle is 30 degrees
        diag1 = power((2.0 * gdf[value_field]) / (power(3, 0.5) / 3.0), 0.5)
        diag2 = diag1 * (power(3, 0.5) / 3.0)
        diag1 = diag1 * scale_factor
        diag2 = diag2 * scale_factor
        xy = [[[cents[i][0] - diag2[i] / 2.0, cents[i][1]],
               [cents[i][0], cents[i][1] + diag1[i] / 2.0], 
               [cents[i][0] + diag2[i] / 2.0, cents[i][1]], 
               [cents[i][0], cents[i][1] - diag1[i] / 2.0]] for i in range(len(diag1))]
        # Create patches
        patches = [Polygon(xy[i]) for i in range(len(xy))]

    elif symbol in ['pentagram']:
        rad = power(gdf[value_field] /
                    ((10.0 * tan(pi / 10.0)) /
                     (3.0 - power(pi / 10.0, 2.0))), 0.5)
        h = 2.0 * ((gdf[value_field] / 10.0) / rad)
        s = h / sin(18.0 * pi / 180.0)
        rad *= scale_factor
        h *= scale_factor
        s *= scale_factor
        # Manually figure out vertex positions.
        xy = []
        for j in range(len(rad)):
            v0 = [cents[j][0], cents[j][1] + rad[j]]
            v1 = [cents[j][0] + h[j], cents[j]
                  [1] + h[j] / tan(36.0 * pi / 180.0)]
            v2 = [cents[j][0] + h[j] + s[j], cents[j]
                  [1] + h[j] / tan(36.0 * pi / 180.0)]
            v3 = [v2[0] - s[j] * sin(54.0 * pi / 180.0),
                  v2[1] - s[j] * cos(54.0 * pi / 180.0)]
            v4 = [cents[j][0] + rad[j] * sin(36.0 * pi / 180.0),
                  cents[j][1] - rad[j] * cos(36.0 * pi / 180.0)]
            v5 = [cents[j][0], cents[j][1] - h[j] / sin(36.0 * pi / 180.0)]
            v6 = [cents[j][0] - rad[j] * sin(36.0 * pi / 180.0),
                  cents[j][1] - rad[j] * cos(36.0 * pi / 180.0)]
            v8 = [cents[j][0] - h[j] - s[j], cents[j]
                  [1] + h[j] / tan(36.0 * pi / 180.0)]
            v7 = [v8[0] + s[j] * sin(54.0 * pi / 180.0),
                  v8[1] - s[j] * cos(54.0 * pi / 180.0)]
            v9 = [cents[j][0] - h[j], cents[j]
                  [1] + h[j] / tan(36.0 * pi / 180.0)]
            xy.append([v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v0])
        patches = [Polygon(xy[i]) for i in range(len(xy))]

    elif symbol in ['hexagram']:
        s = power(gdf[value_field] / 3.0 * power(3, 0.5), 0.5)
        h = 2.0 * power(0.75 * power(s, 2), 0.5)
        s *= scale_factor
        h *= scale_factor
        # Manually figure out vertex positions
        xy = []
        for j in range(len(s)):
            v0 = [cents[j][0], cents[j][1] + h[j]]
            v1 = [cents[j][0] + s[j] / 2.0, cents[j][1] + h[j] / 2.0]
            v2 = [cents[j][0] + 1.5 * s[j], cents[j][1] + h[j] / 2.0]
            v3 = [cents[j][0] + s[j], cents[j][1]]
            v4 = [cents[j][0] + 1.5 * s[j], cents[j][1] - h[j] / 2.0]
            v5 = [cents[j][0] + s[j] / 2.0, cents[j][1] - h[j] / 2.0]
            v6 = [cents[j][0], cents[j][1] - h[j]]
            v7 = [cents[j][0] - s[j] / 2.0, cents[j][1] - h[j] / 2.0]
            v8 = [cents[j][0] - 1.5 * s[j], cents[j][1] - h[j] / 2.0]
            v9 = [cents[j][0] - s[j], cents[j][1]]
            v10 = [cents[j][0] - 1.5 * s[j], cents[j][1] + h[j] / 2.0]
            v11 = [cents[j][0] - s[j] / 2.0, cents[j][1] + h[j] / 2.0]
            xy.append([v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v0])
        patches = [Polygon(xy[i]) for i in range(len(xy))]

    else:
        print("symbol parameter invalid, using 'circle', use help(propsymbols) for all options")
        # Work out circular symbol radii
        radii = (power(gdf[value_field] / pi, 0.5) * scale_factor).tolist()
        # Create patches based on radii
        patches = [Circle(cents[i], radius=radii[i])
                   for i in range(len(radii))]

    return PatchCollection(
        patches, zorder=10, edgecolors=edgecolors, facecolors=facecolors, alpha=alpha)
