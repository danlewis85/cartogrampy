from __future__ import print_function, division

def dcn_cartogram(geodataframe,value_field,id_field=None,geom_field='geometry',iterations=99,verbose=True):
    
    """Returns a Dougenik, Chrisman, Niemayer cartogram for a given GeoPandas GeoDataFrame as a new GeoDataFrame.

    Args:
        geodataframe (geopandas.geodataframe.GeoDataFrame): Input GeoDataFrame.
        value_field: Field name of values to be used to produce DCN cartogram.
        id_field (str, optional): Field name of values uniquely identifying rows in geodataframe, defaults to None.
        geom_field (str, optional): Field name of geometry column in input GeoDataFrame, defaults to 'geometry'.
        iterations (int, optional): Number of iteration fo the algorithm to perform. Default of 99 is arbitrary.
        verbose (bool, optional): Prints the iteration and mean and max error, default = True.
        
    Returns:
        GeoPandas GeoDataFrame.
    """

    # imports
    from numpy import mean, pi, power, maximum, minimum, array, concatenate, lexsort, any, str, argsort, searchsorted, split, take
    from numpy.linalg import norm
    from shapely.geometry import Polygon
    from geopandas import GeoDataFrame

    if id_field:
        gdf = geodataframe[[geom_field,value_field,id_field]].copy()
    else:
        gdf = geodataframe[[value_field,geom_field]].copy()
        gdf['id_field'] = gdf.index
        id_field = 'id_field'
    
    # compute sum of value_field and store
    TotalValue = gdf[value_field].sum()
    # Main loop for iterations.
    for i in range(iterations):
        # Calculate all the current statistics for gdf
        # Total area of all geometries for current iteration
        TotalArea = gdf[geom_field].area.sum()
        # Desired areas relative to total area of current iteration.
        Desired = TotalArea * (gdf[value_field]/TotalValue)
        # Circular radius for each geometry as a function of area of current iteration
        Radius = power(gdf[geom_field].area/pi,0.5)
        # Mass as difference in desired and actual radii 
        Mass = power(Desired/pi,0.5) - Radius
        # Ratio of desired to current area as an error rate.
        SizeError = maximum(gdf[geom_field].area,Desired)/minimum(gdf[geom_field].area,Desired)
        # Work out force reduction factor based on mean error.
        ForceReductionFactor = 1.0 / (1.0 + mean(SizeError))
        # Create centroids for convenience. Need to be multi-polygon centroids if multipolygon.
        centroids = gdf[geom_field].centroid
        
        # Now prepare the geodataframes - ensure singlepart polygons and deduplicate geometry points.

        # If multipolygons exist (Multipart), convert to polygons (Singlepart).
        multi = False
        if gdf[gdf[geom_field].type == 'MultiPolygon'][geom_field].any():
            # Convert to geodataframe with singlepolys.
            # First make dictionary to hold data.
            data = {}
            for col in gdf.columns:
                if col != geom_field:
                    data[col] = []
            
            # Now iterate over gdf to get data and geometry rows.
            geom = []
            for idx,row in gdf.iterrows():
                if row[geom_field].type == 'Polygon':
                    for col in data.keys():
                        data[col].append(row[col])
                    geom.append(row[geom_field])
                elif row[geom_field].type == 'MultiPolygon':
                    for g in row[geom_field]:
                        for col in data.keys():
                            data[col].append(row[col])
                        geom.append(g)
                else:
                    print("Geometry Error")
            
            # Reset gdf to the singlepart feature geodataframe
            gdf = GeoDataFrame(data,geometry = geom)        
            # Flag that the geometry in question is originally multipart.
            multi = True
        
        # Now that we have singlepart features, deduplicate geometry points.

        # creates a series of numpy arrays of boundary points
        pnts = gdf[geom_field].apply(lambda x: array(x.exterior))

        # create a lookup based on count of points per polygon
        cnts = pnts.apply(lambda x: len(x))

        # Concatenate all points into a single long list of points.
        pnts = concatenate(array(pnts))

        # Calculate unique points in full dataset.
        # From: http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
        ind = lexsort(pnts.T)
        upnts = pnts[ind[concatenate(([True],any(pnts[ind[1:]]!=pnts[ind[:-1]],axis=1)))]]

        # This creates an id lookup for pnts to upnts.
        # It creates string views of the data, which is WAY quicker to deal with that 2d arrays.
        # Using http://stackoverflow.com/questions/32191029/getting-the-indices-of-several-elements-in-a-numpy-array-at-once
        a = pnts.ravel().view((str, pnts.itemsize * pnts.shape[1]))
        b = upnts.ravel().view((str, upnts.itemsize * upnts.shape[1]))
        sorter = argsort(b)
        idv = sorter[searchsorted(b, a, sorter=sorter)]
        
        # This is the main algorithm, and the bit that really takes the time.
        points = []
        for xy in upnts:
            x = x0 = xy[0]
            y = y0 = xy[1]
            idx = 0
            for cxy in centroids.apply(lambda x: [x.x,x.y]):
                dist = norm(array(cxy)-array([x,y]))
                if dist > Radius[idx]:
                    Fij = Mass[idx] * Radius[idx] / dist
                else:
                    xF = dist / Radius[idx]
                    Fij = Mass[idx] * (power(xF,2)) * (4 - (3 * xF))
                Fij = Fij * ForceReductionFactor / dist
                
                x = (x0 - cxy[0]) * Fij + x
                y = (y0 - cxy[1]) * Fij + y
                idx +=1
            points.append([x,y])
        points = array(points)
        # Reconstruct the full points list from the unique points using take and the index list.
        # Remake the array - split by counts lookup.
        # NB without the -1 index you get a 0 length array at the end.
        repnts = split(take(points,idv,axis=0),cnts.cumsum())[:-1]
        new_geom = [Polygon(p) for p in repnts]
        gdf = GeoDataFrame(gdf[[value_field,id_field]],geometry=new_geom)
        
        if multi:
            # Dissolve back to original multipart polygon.
            gdf = gdf.dissolve(by = id_field,as_index=False)
        if verbose:
            print("iteration:", i+1, "Mean Error:", mean(maximum(gdf[geom_field].area,Desired)/minimum(gdf[geom_field].area,Desired)),"Max Error:",max(maximum(gdf[geom_field].area,Desired)/minimum(gdf[geom_field].area,Desired)))
    return gdf