import geopandas as gpd
import pandas as pd
import matplotlib
matplotlib.use("qt5agg")
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from cartogrampy.components.cartogram import Cartogram

# ldn_boro = gpd.read_file('data/LDN_Boro.geojson').to_crs(epsg=27700)
# ldn_pop  = pd.read_excel('data/Pandas_Lon_Pop.xlsx')
# data = ldn_boro.merge(ldn_pop, how='left', left_on='GSS_CODE', right_on='New Code')
# c = Cartogram(data, value_field=2015 )
# dcn = c.dcn()
# dorling = c.dorling(iterations=200,verbose=True)
# ncc = c.noncont()

# Get some data for US States
usstate = gpd.read_file('data/US_State_2016_5m.geojson')
# Set crs
usstate.crs = {'init': u'epsg:4269'}
# Get continental US and project to NAD83 Contiguous US Albers.
usstate = usstate[~usstate['STUSPS'].isin(['AK', 'HI', 'AS', 'PR',
                                           'GU', 'MP', 'VI'])].to_crs({'init': 'epsg:5070'})
# Read in state populations
state_pop = pd.read_excel('data/Pandas_US_Pop.xlsx')

# Merge population data
usstate = usstate.merge(state_pop, how='left', left_on='STUSPS', right_on='Code')

c = Cartogram(usstate, value_field=2016)
dcn = c.dcn(30)

# Map the outcome of dcn
# fig, ax = plt.subplots(1,2,figsize=(16,8))
# usstate.plot(column = 2016,cmap = 'PuRd', legend=True, ax=ax[0])
# dcn.plot(color='r',ax=ax[1])
# ax[0].axis('equal')
# ax[1].axis('equal')
# plt.show()

# Map the outcome of dorling
# fig, ax = plt.subplots(1,2,figsize=(16,8))
# data.plot(column = 2015, cmap = 'PuRd', scheme = 'Quantiles', k=4,legend=True, ax=ax[0])
# dorling.plot(color='r',ax=ax[1])
# ax[0].axis('equal')
# ax[1].axis('equal')
# plt.show()

# Map the outcome of NCC
# fig, ax = plt.subplots(1,figsize=(8,8))
# data.plot(color='black',ax=ax,alpha = 0.8,zorder=0)
# ncc.plot(color='r',ax=ax)
# ax.axis('equal')
# ax.axis('equal')
# plt.show()
