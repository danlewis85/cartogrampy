import geopandas as gpd
import pandas as pd
# from cartogrampy.components.cartogram import Cartogram

ldn_boro = gpd.read_file('data/LDN_Boro.geojson').to_crs(epsg=27700)
ldn_pop  = pd.read_excel('data/Pandas_Lon_Pop.xlsx')
data = ldn_boro.merge(ldn_pop, how='left', left_on='GSS_CODE', right_on='New Code')
c = Cartogram(data, value_field=2015 )
# c.dcn(verbose=True)
