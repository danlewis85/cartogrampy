from __future__ import division, print_function
# Copyright. 2017. Daniel Lewis.
# https://github.com/danlewis85/cartogrampy/blob/master/LICENSE

def borders_from_dataframe(df, idVariable=None,geom_field = 'geometry'):
        
    """Returns a PySAL weights object in which weights are lengths of shared border.

    Args:
        df (geopandas.geodataframe.GeoDataFrame): Input GeoDataFrame.
        idVariance (optional): Field name of values uniquely identifying rows in geodataframe, defaults to None.
        geom_field (str, optional): Field name of geometry column in input GeoDataFrame, defaults to 'geometry'.
    
    Returns:
        PySAL W weights object.
    """
        
    # Function to compute neighbors and length of shared borders
    # Based on: https://github.com/volaya/processing_pysal/blob/master/ext-libs/pysal/contrib/shared_perimeter_weights.py
        
    # imports
    from pysal.weights.Contiguity import Rook, W
    from shapely.geometry import asShape

    rook = Rook.from_dataframe(df,idVariable=idVariable)
    polygons = df.set_index(idVariable)[geom_field].apply(lambda x: asShape(x)).to_dict()
    new_weights = {}
    for i in rook.neighbors:
        a = polygons[i]
        new_weights[i] = [a.intersection(polygons[j]).length for j in rook.neighbors[i]]
    return W(rook.neighbors,new_weights)

def dorling_cartogram(geodataframe, value_field, id_field = None, geom_field = 'geometry',
                     ratio = 0.4,friction = 0.25,iterations = 99,verbose=True):
    
    """Returns a Dorling cartogram for a given GeoPandas GeoDataFrame as a new GeoDataFrame.

    Args:
        geodataframe (geopandas.geodataframe.GeoDataFrame): Input GeoDataFrame.
        value_field: Field name of values to be used to produce Dorling cartogram.
        id_field (str, optional): Field name of values uniquely identifying rows in geodataframe, defaults to None.
        geom_field (str, optional): Field name of geometry column in input GeoDataFrame, defaults to 'geometry'.
        ratio (float, optional): Ratio by which to combine repulsion and attraction forces, default = 0.4
        friction (float, optional): Relaxation factor for forces, default = 0.25
        iterations (int, optional): Number of iteration fo the algorithm to perform. Default of 99 is arbitrary.
        verbose (bool, optional): Prints the iteration and displacement, default = True.
        
    Returns:
        GeoPandas GeoDataFrame.
    """
    from geopandas import GeoDataFrame
    from pandas import DataFrame
    from numpy import sum, pi, power
    from shapely.affinity import translate
    
    # Copy the required data from the input geodataframe.
    if id_field:
        df = geodataframe[[id_field,value_field,geom_field]].copy()
    else:
        df = geodataframe[[value_field,geom_field]].copy()
        df['id_field'] = df.index
        id_field = 'id_field'
    
    # Get lengths of shared borders.
    Wp = borders_from_dataframe(df, id_field, geom_field)
    
    # Get dictionary of perimeters
    perimeter = df.set_index(id_field)[geom_field].length.to_dict()
    
    # Now convert df to geodataframe of centroids
    df = GeoDataFrame(df,geometry=df[geom_field].centroid)
    
    # Work out scale and radii - seems inefficient to deal with the geometries in this way.
    # No idea how to write this is a PEP8 compliant way.
    total_dist = sum([df.loc[df[id_field] == i,geom_field].centroid.tolist()[0].distance(df.loc[df[id_field] == j,geom_field].centroid.tolist()[0]) for i in Wp.neighbors for j in Wp[i]])
    total_radius = sum([power(df.loc[df[id_field] == i,value_field].values[0]/pi,0.5) + power(df.loc[df[id_field] == j,value_field].values[0]/pi,0.5) for i in Wp.neighbors for j in Wp[i]])

    scale = total_dist / total_radius
    
    radius = df.set_index(id_field)[value_field].apply(lambda x: scale * power(x/pi,0.5)).to_dict()
    widest = radius[max(radius, key=lambda i: radius[i])]
    
    if verbose:
        print('Scaling by', scale, 'widest is', widest)
    
    # Start the main loop
    for i in range(iterations):
        displacement = 0.0
        for idx, row in df.iterrows():
            # Set up for values for each iteration.
            xrepel = 0.0
            yrepel = 0.0
            xattract = 0.0
            yattract = 0.0
            closest = widest
        
            # Get neighboring circles.
            distband = row[geom_field].centroid.buffer(widest + radius[row[id_field]])
            neighbors = df[df[geom_field].centroid.within(distband)]
            # Remove self from neighbors dataframe
            neighbors = neighbors[neighbors[id_field] != row[id_field]]
            
            # For neighbors calculate the attractive and repulsive forces acting.
            if len(neighbors) > 0:
                
                for nidx, nrow in neighbors.iterrows():
                    dist = row[geom_field].centroid.distance(nrow[geom_field].centroid)
                    if dist < closest:
                        closest = dist
                    overlap = radius[row[id_field]] + radius[nrow[id_field]] - dist
                    # Calculate repelling forces for intersecting circles
                    if overlap > 0.0:
                        if dist > 1.0:
                            xrepel -= overlap * (nrow[geom_field].centroid.x - row[geom_field].centroid.x) / dist
                            yrepel -= overlap * (nrow[geom_field].centroid.y - row[geom_field].centroid.y) / dist
                    # Calculate attractive forces for circles
                    if overlap < 0.0:
                        try:
                            overlap = abs(overlap) * Wp[row[id_field]][nrow[id_field]]/perimeter[row[id_field]]
                        except KeyError:
                            gap = 0.0
                        xattract = xattract + overlap * (nrow[geom_field].centroid.x - row[geom_field].centroid.x) / dist
                        yattract = yattract + overlap * (nrow[geom_field].centroid.y - row[geom_field].centroid.y) / dist

            # Calculate combined effect of attraction and repulsion
            atrdst = power(power(xattract, 2) + power(yattract, 2), 0.5)
            repdst = power(power(xrepel, 2) + power(yrepel, 2), 0.5)
            if repdst > closest:
                xrepel = closest * xrepel / (repdst + 1.0)
                yrepel = closest * yrepel / (repdst + 1.0)
                repdst = closest
            if repdst > 0:
                xtotal = (1.0 - ratio) * xrepel + ratio * (repdst * xattract / (atrdst + 1.0))
                ytotal = (1.0 - ratio) * yrepel + ratio * (repdst * yattract / (atrdst + 1.0))
            else:
                if atrdst > closest:
                    xattract = closest * xattract / (atrdst + 1.0)
                    yattract = closest * yattract / (atrdst + 1.0)
                xtotal = xattract
                ytotal = yattract
            displacement += power(power(xtotal, 2) + power(ytotal, 2), 0.5)
            
            # Record the vectors
            xvector = friction * xtotal
            yvector = friction * ytotal
            
            # update position
            df.loc[idx,geom_field] = translate(row[geom_field], xoff = xvector, yoff = yvector)
        
        displacement = displacement/len(df)

        if verbose:
            print("iter:",i,"displacement:", displacement)
    
    df = df.merge(DataFrame(radius.items(), columns=[id_field, 'Radius']),on = id_field)
    return GeoDataFrame(df,geometry=[df.loc[b,'geometry'].buffer(df.loc[b,'Radius']) for b in range(len(df))])