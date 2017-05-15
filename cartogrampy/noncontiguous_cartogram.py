def noncontiguous_cartogram(geodataframe, value_field, id_field = None, geom_field = 'geometry', 
                            position = 'centroid', anchor_rank = 1, anchor_id = None):
    
    """Returns a non-contiguous cartogram for a given GeoPandas GeoDataFrame as a new GeoDataFrame

    Args:
        geodataframe (geopandas.geodataframe.GeoDataFrame): Input GeoDataFrame.
        value_field (str): Field name of values to be used to produce non-contiguous cartogram.
        id_field (str, optional): Field name of values uniquely identifying rowsin geodataframe, defaults to None.
        geom_field (str, optional): Field name of geometry column in input GeoDataFrame, defaults to 'geometry'.
        position (str, optional): Center position of scaled polygons, defaults to centroid.
            'centroid' - GeoDataFrame centroid property, shape centroid.
            'center' - Centroid of GeoDataFrame envelope (bounding box) property, envelope centroid.
            'representative point' - Point returned by GeoDataFrame represntative point method, guaranteed to be inside polygon.
        anchor_rank (int, optional): Rank of row value-density to be used as reference (anchor) for scaling of geometry. Default = 1.
        anchor_id (optional): Value of row to use as anchor, 'id_field' must also be specified, defaults to None.

    Returns:
        GeoPandas GeoDataFrame
    """
    from numpy import power
    from shapely.affinity import scale
    from geopandas import GeoDataFrame
    
    # Make a copy of the geodataframe containing only the relevant fields.
    if id_field:
        gdf = geodataframe[[value_field,id_field,geom_field]].copy()
    else:
        gdf = geodataframe[[value_field,geom_field]].copy()
    
    # Calculate geometry positions based on
    if position.lower() in ['centroid']:
        gdf['cent'] = gdf[geom_field].centroid
    elif position.lower() in ['center', 'centre']:
        gdf['cent'] = gdf[geom_field].envelope.centroid
    elif position.lower() in ['representative point', 'rep']:
        gdf['cent'] = gdf[geom_field].representative_point()
    else:
        # If position parameter not recognised, default to centroid.
        print("position parameter invalid, using 'centroid'. Options are 'centroid','center', or 'representative point'")
        gdf['cent'] = gdf[geom_field].centroid
    
    # Work out value densities and ranks.
    gdf['density'] = gdf[value_field]/gdf.area
    gdf['rank'] = gdf['density'].rank(axis=0,method='first',ascending=False)
    
    # Get the appropriate anchor depending on whether anchor_id or anchor_rank have been specified
    if anchor_id:
        if id_field:
            if anchor_id in gdf[id_field].values:
                anchor = gdf[gdf[id_field] == anchor_id]['density'].values[0]
            else:
                print("anchor_id not recognized in id_field, defaulting to anchor_rank = 1")
                anchor = gdf[gdf['rank'] == 1]['density'].values[0]
        else:
            print("id_field not specified for anchor_id, defaulting to anchor_rank = 1")
            anchor = gdf[gdf['rank'] == 1]['density'].values[0]
    else:
        anchor = gdf[gdf['rank'] == anchor_rank]['density'].values[0]
    
    # Work out the scaling for each polygon
    gdf['scale'] = (1.0/power(anchor,0.5)) * power(gdf[value_field]/gdf.area,0.5) 
    
    # NB affine transformations are linear
    new_geom = [scale(g[1][geom_field],xfact=g[1]['scale'],yfact=g[1]['scale'],origin=g[1]['cent']) for g in gdf.iterrows()]
    
    # Clean up
    del gdf['density']
    del gdf['rank']
    del gdf['cent']
    
    return GeoDataFrame(gdf,geometry=new_geom)