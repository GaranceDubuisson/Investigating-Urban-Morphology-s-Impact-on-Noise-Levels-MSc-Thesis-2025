#%% import

import osmnx as ox
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import box
from shapely.geometry import mapping
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

import os
from openpyxl import Workbook




#%% tags

tags_green = {'leisure': ['park', 'garden'], 'landuse': ['forest', 'grass', 'cemetery'], 'natural': 'wood'}
tags_water = {
    'natural': ['water', 'coastline', 'sea'],  # Lacs, mers, côtes
    'waterway': ['river', 'stream', 'canal'],  # Rivières, ruisseaux, canaux
    'landuse': ['reservoir', 'basin'],  # Réservoirs et bassins artificiels
    'leisure': 'marina'  # Facultatif : ports de plaisance
}
tags_buildings = {'building': True}

tags_public = {
    'building': [
        'hospital',    # Hôpitaux
        'school',      # Écoles
        'university',  # Universités
        'public',      # Bâtiments publics
        'civic',       # Bâtiments civiques
        #'office',      # Bureaux administratifs
        'yes'          # Autres bâtiments généraux
    ]
}

tags_office = {
    'building': ['office']
}

tags_commercial = {
        'shop': True,
        'amenity': ['restaurant', 'cafe', 'bar', 'fast_food', 'bank', 'post_office'],
        'office': True,
        'building': ['retail', 'commercial', 'supermarket'],
        'landuse': 'commercial'
    }


#%% function intersections count

def count_streets_per_node_manual(G):
    # Dictionnaire pour stocker le nombre de rues par nœud
    streets_per_node = {}
    
    for node, data in G.nodes(data=True):
        # Récupérer le nombre d'arêtes (routes) connectées à ce nœud
        num_streets = len(list(G.edges(node)))
        streets_per_node[node] = num_streets
    
    return streets_per_node


#%% function mean compactness

def calculate_mean_compactness(buildings):
    buildings = buildings[buildings.is_valid]
    
    # Calculer la surface et le périmètre pour chaque bâtiment
    buildings['area'] = buildings.geometry.area
    buildings['perimeter'] = buildings.geometry.length

    buildings = buildings[(buildings['area'] > 0) & (buildings['perimeter'] > 0)]
    # Calculer le compactness ratio pour chaque bâtiment
    buildings['compactness'] = (4 * np.pi * buildings['area']) / (buildings['perimeter'] ** 2)
    
    # Calculer la moyenne du ratio de compacité, en excluant les valeurs NaN
    # mean_compactness = buildings['compactness'].mean()
    
    valid_buildings = buildings[buildings['perimeter'] > 0]

    # Si aucun bâtiment valide, retourner NaN
    if valid_buildings.empty:
        return np.nan
    
    # Calculez la moyenne pondérée
    mean_compactness = (valid_buildings['compactness'] * valid_buildings['area']).sum() / valid_buildings['area'].sum()
    
    return mean_compactness

#%% SW/H

def calculate_sw_h(buildings: gpd.GeoDataFrame, streets: gpd.GeoDataFrame, buffer_distance: float, mean_height: float) -> gpd.GeoDataFrame:
    
    streets['buffer'] = streets.geometry.buffer(buffer_distance)

    streets['left_building_distance'] = streets['buffer'].apply(
        lambda buffer: buildings[buildings.geometry.intersects(buffer)].distance(buffer).min()
        if not buildings[buildings.geometry.intersects(buffer)].empty else np.nan)
    streets['street_width'] = 2 * streets['left_building_distance']
    streets['street_width'] = streets['street_width'].fillna(10)
    streets['street_width'] = streets['street_width'].replace(0, 10)

    streets['mean_building_height'] = streets['buffer'].apply(
        lambda buffer: buildings[buildings.geometry.intersects(buffer)]['height'].mean()
        if not buildings[buildings.geometry.intersects(buffer)].empty else np.nan)
    streets['mean_building_height'].fillna(mean_height, inplace=True)

    streets['SW_H'] = np.where(
        (streets['mean_building_height'] > 0) & (streets['street_width'] > 0),
        streets['street_width'] / streets['mean_building_height'],
        np.nan)

    if streets['SW_H'].notna().any():
        global_sw_h = np.average(streets['SW_H'].dropna(), weights=streets.loc[streets['SW_H'].notna(), 'geometry'].length)
        return global_sw_h
    else:
        return np.nan

#%% function for exceptions

def get_features_with_exception_from_polygon(polygon, tags):
    """
    Extrait les entités OSM correspondant aux tags spécifiés à partir d'un polygone,
    en gérant les exceptions éventuelles.
    
    Paramètres :
    - polygon : shapely.geometry.polygon.Polygon
        Le polygone à partir duquel extraire les entités.
    - tags : dict
        Un dictionnaire de tags OSM pour filtrer les entités (par ex., {'building': True}).
        
    Retour :
    - GeoDataFrame contenant les entités correspondant aux tags, ou None si une exception se produit.
    """
    try:
        features = ox.features_from_polygon(polygon, tags)
        gdf = gpd.GeoDataFrame(features, crs='EPSG:4326').to_crs(epsg=3857)
        
        if gdf.empty:
            print("Aucune entité trouvée pour les tags donnés.")
            return None
        
        return gdf
    
    except Exception as e:
        # En cas d'erreur, afficher un message et renvoyer None
        print(f"Erreur lors de l'extraction des entités pour les tags {tags}: {e}")
        return None

def get_features_with_exception_from_point(coord, tags, dist):
    """
    Extrait les entités OSM correspondant aux tags spécifiés à partir d'un polygone,
    en gérant les exceptions éventuelles.
    
    Paramètres :
    - polygon : shapely.geometry.polygon.Polygon
        Le polygone à partir duquel extraire les entités.
    - tags : dict
        Un dictionnaire de tags OSM pour filtrer les entités (par ex., {'building': True}).
        
    Retour :
    - GeoDataFrame contenant les entités correspondant aux tags, ou None si une exception se produit.
    """
    try:
        features = ox.features_from_point(coord, tags, dist)
        gdf = gpd.GeoDataFrame(features, crs='EPSG:4326').to_crs(epsg=3857)
        
        if gdf.empty:
            print("Aucune entité trouvée pour les tags donnés.")
            return None
        
        return gdf
    
    except Exception as e:
        # En cas d'erreur, afficher un message et renvoyer None
        print(f"Erreur lors de l'extraction des entités pour les tags {tags}: {e}")
        return None

#%% For polygons

CPH_path = r'C:\Users\Garance Dubuisson\OneDrive - Danmarks Tekniske Universitet\Thesis\QGIS\v4\CPH_districts.shp'

def calculate_urban_parameters_from_shapefile2(shapefile_path = CPH_path):
    polygons = gpd.read_file(shapefile_path).to_crs(epsg=4326)
    
    results = []
    
    for idx, poly in polygons.iterrows(): 
        polygon = poly.geometry
        #if polygon.is_valid:
            #print("Polygon is valid.")
          
        gmc = poly.get('GMC', None)
        print(f"Processing district with GMC: {gmc}")
        
        # in CRS 3857
        buildings = get_features_with_exception_from_polygon(polygon, tags_buildings)
        water = get_features_with_exception_from_polygon(polygon, tags_water)
        green = get_features_with_exception_from_polygon(polygon, tags_green)
        public = get_features_with_exception_from_polygon(polygon, tags_public)
        office = get_features_with_exception_from_polygon(polygon, tags_office)
        commercial = get_features_with_exception_from_polygon(polygon, tags_commercial)
                
        
        # Calculer les paramètres
        #total_area = polygon.area
        total_area = poly.get('Area_km2')*1000000 #in m2
        tot_area = polygon.area
        pop_dens = poly.get('Pop_densit')
        building_area = buildings.geometry.area.sum()
        GSI = building_area / total_area
        print("GSI from QGIS:", GSI, "vs with Python:", building_area/tot_area)
        buildings['building_area'] = buildings.geometry.area
        if 'building:levels' in buildings.columns:
            buildings['building:levels'] = pd.to_numeric(buildings['building:levels'], errors='coerce')
            av_levels = buildings['building:levels'].mean()
            #print('av levels:', av_levels)
            buildings['levels'] = buildings['building:levels'].fillna(av_levels).astype(float)
        else:
            buildings['levels'] = 1.0
        buildings['floor_area'] = buildings['building_area'] * buildings['levels']
        floor_area = buildings['floor_area'].sum()
        FSI = floor_area / total_area
        porosity = 1 - GSI
        print("porosity with QGIS:", porosity, "vs with Python:", 1-building_area/tot_area)
        commercial_area = commercial.geometry.area.sum() if (not commercial is None) else 0
        Com_ratio = commercial_area / building_area if building_area > 0 else 0
        public_area = public.geometry.area.sum() if (not public is None) else 0
        Pub_ratio = public_area / building_area if building_area > 0 else 0
        office_area = office.geometry.area.sum() if (not office is None) else 0
        #office_area = office.area.sum() if (not office is None) else 0
        Off_ratio = office_area / building_area if building_area > 0 else 0
        green_area = green.geometry.area.sum() if (not green is None) else 0
        #green_area = green.area.sum() if (not green is None) else 0
        green_ratio = green_area / total_area
        water_area = water.geometry.area.sum() if (not water is None) else 0
        water_ratio = water_area / total_area

        # Environnement bâti
        buildings['height'] = buildings['levels'] * 3.2
        mean_height = np.average(buildings['height'], weights=buildings['building_area'])
        mean_compactness = calculate_mean_compactness(buildings)
        num_buildings = len(buildings)

        # Réseau de routes
        streets = ox.graph_from_polygon(polygon, network_type='drive')
        intersections = count_streets_per_node_manual(streets)
        num_intersections = len([node for node, count in intersections.items() if count > 1])
        #print('intersections:', num_intersections)
        streets = ox.graph_to_gdfs(streets, nodes=False, edges=True).to_crs(epsg=3857)
        buffer = 20
        SWH = calculate_sw_h(buildings, streets, buffer, mean_height)
        road_length = streets.geometry.length.sum()
        road_coverage = road_length / total_area
        streets['road_area'] = streets['street_width'] * streets['length']
        road_coverage_best = streets['road_area'].sum() / total_area

        
        # Also add the acoustics and traffic
        L_den = poly.get('L_den')
        L_night = poly.get('L_night')
        AADT = poly.get('AADT')
        if AADT is not None:
            AADT=AADT/total_area # veh/day/m²          
        Total_AADT = poly.get('Total_AADT')
        
        # Stocker les résultats
        results.append({
            'District id': gmc,
            'Total area [km2]': total_area/1000000,
            'Population density [p/km2]': pop_dens,
            'GSI (Coverage)': GSI,
            'FSI (Intensity)': FSI,
            'Porosity': porosity,
            'Road Coverage': road_coverage_best,
            'Green Space Ratio': green_ratio,
            'Water Surface Ratio': water_ratio,
            'Commercial Use Ratio': Com_ratio,
            'Public Use Ratio': Pub_ratio,
            'Office Use Ratio': Off_ratio,
            'Mean compactness': mean_compactness,
            'Mean Height [m]': mean_height,
            'Sky View Factor SW/H': SWH,
            'Number of Buildings [/km2]': num_buildings*1000000/total_area,
            'Road Length [m/km2]': road_length*1000000/total_area,
            'Number of Intersections [/km2]': num_intersections*1000000/total_area,
            'Total AADT [veh/day/km2]': Total_AADT*1000000/total_area,
            'L_den [dB]': L_den,
            'L_night [dB]': L_night
        })
    
    return results

#%% Modif for better area

CPH_path = r"C:\Users\Garance Dubuisson\OneDrive - Danmarks Tekniske Universitet\Thesis\QGIS\v4\CPH_districts2.shp"


def calculate_urban_parameters_from_shapefile(shapefile_path = CPH_path):
    polygons = gpd.read_file(shapefile_path).to_crs(epsg=3857)
    polygons["tot_area"] = polygons.geometry.area
    
    results = []
    polygons = polygons.to_crs(epsg=4326)
    
    for idx, poly in polygons.iterrows(): 
        polygon = poly.geometry
        #if polygon.is_valid:
            #print("Polygon is valid.")
          
        gmc = poly.get('GMC', None)
        print(f"Processing district with GMC: {gmc}")
        
        # in CRS 3857
        buildings = get_features_with_exception_from_polygon(polygon, tags_buildings)
        water = get_features_with_exception_from_polygon(polygon, tags_water)
        green = get_features_with_exception_from_polygon(polygon, tags_green)
        public = get_features_with_exception_from_polygon(polygon, tags_public)
        office = get_features_with_exception_from_polygon(polygon, tags_office)
        commercial = get_features_with_exception_from_polygon(polygon, tags_commercial)
                
        
        # Calculer les paramètres
        total_area = poly.get('tot_area')
        pop_dens = poly.get('Pop_dens')
        if gmc == 9 or gmc == 10: # wrong data in Osterbro
            pop_dens = None
        building_area = buildings.geometry.area.sum()
        GSI = building_area / total_area
        buildings['building_area'] = buildings.geometry.area
        if 'building:levels' in buildings.columns:
            buildings['building:levels'] = pd.to_numeric(buildings['building:levels'], errors='coerce')
            av_levels = buildings['building:levels'].mean()
            #print('av levels:', av_levels)
            buildings['levels'] = buildings['building:levels'].fillna(av_levels).astype(float)
        else:
            buildings['levels'] = 1.0
        buildings['floor_area'] = buildings['building_area'] * buildings['levels']
        floor_area = buildings['floor_area'].sum()
        FSI = floor_area / total_area
        porosity = 1 - GSI
        commercial_area = commercial.geometry.area.sum() if (not commercial is None) else 0
        Com_ratio = commercial_area / building_area if building_area > 0 else 0
        public_area = public.geometry.area.sum() if (not public is None) else 0
        Pub_ratio = public_area / building_area if building_area > 0 else 0
        office_area = office.geometry.area.sum() if (not office is None) else 0
        #office_area = office.area.sum() if (not office is None) else 0
        Off_ratio = office_area / building_area if building_area > 0 else 0
        green_area = green.geometry.area.sum() if (not green is None) else 0
        #green_area = green.area.sum() if (not green is None) else 0
        green_ratio = green_area / total_area
        water_area = water.geometry.area.sum() if (not water is None) else 0
        water_ratio = water_area / total_area

        # Environnement bâti
        buildings['height'] = buildings['levels'] * 3.2
        mean_height = np.average(buildings['height'], weights=buildings['building_area'])
        mean_compactness = calculate_mean_compactness(buildings)
        num_buildings = len(buildings)

        # Réseau de routes
        streets = ox.graph_from_polygon(polygon, network_type='drive')
        intersections = count_streets_per_node_manual(streets)
        num_intersections = (len([node for node, count in intersections.items() if count > 1])*1000000)/total_area
        #print('intersections:', num_intersections)
        streets = ox.graph_to_gdfs(streets, nodes=False, edges=True).to_crs(epsg=3857)
        buffer = 20
        SWH = calculate_sw_h(buildings, streets, buffer, mean_height)
        road_length = streets.geometry.length.sum()
        road_coverage = road_length / total_area
        streets['road_area'] = streets['street_width'] * streets['length']
        road_coverage_best = streets['road_area'].sum() / total_area

        
        # Also add the acoustics and traffic
        L_den = poly.get('L_den')
        L_night = poly.get('L_night')
        Total_AADT = poly.get('Total_AADT')
        if Total_AADT is not None:
            AADT=Total_AADT*1000000/total_area # veh/day/km²
        if Total_AADT == 0:
            AADT=None
        
        # Stocker les résultats
        results.append({
            'District id': gmc,
            'Total area [km2]': total_area/1000000,
            'Population density [p/km2]': pop_dens,
            'GSI (Coverage)': GSI,
            'FSI (Intensity)': FSI,
            'Porosity': porosity,
            'Road Coverage [m/m2]': road_coverage,
            'Road Coverage [n.u.]': road_coverage_best,
            'Green Space Ratio': green_ratio,
            'Water Surface Ratio': water_ratio,
            'Commercial Use Ratio': Com_ratio,
            'Public Use Ratio': Pub_ratio,
            'Office Use Ratio': Off_ratio,
            'Mean compactness': mean_compactness,
            'Mean Height [m]': mean_height,
            'Sky View Factor SW/H': SWH,
            'Number of Buildings': num_buildings,
            'Road Length [m]': road_length,
            'Number of Intersections [/km2]': num_intersections,
            'AADT [veh/day/km2]': AADT,
            'L_den [dB]': L_den,
            'L_night [dB]': L_night
        })
    
    return results

#%% Pattern

def get_pattern_polygons(shapefile_path, coord, dist=5000):
    # Charger le shapefile et les données OSM
    polygons = gpd.read_file(shapefile_path).to_crs(epsg=4326)
    polygons = polygons.to_crs(epsg=3857)  # Reprojection en EPSG:3857 pour l'affichage

    # Charger les couches de rues, bâtiments, eau et espaces verts
    streets = ox.graph_from_point(coord, dist=dist, dist_type='bbox', network_type='drive', simplify=True)
    streets = ox.project_graph(streets)  # Reprojection du graphe en EPSG:3857

    buildings = get_features_with_exception_from_point(coord, tags=tags_buildings, dist=dist)
    water = get_features_with_exception_from_point(coord, tags=tags_water, dist=dist)
    green = get_features_with_exception_from_point(coord, tags=tags_green, dist=dist)


    # # Graphique 1 : Bâtiments, cours d'eau, espaces verts, avec polygones
    fig1, ax1 = plt.subplots(figsize=(10, 10), facecolor='black')
    buildings.plot(ax=ax1, facecolor='black', edgecolor='black')
    # green.plot(ax=ax1, facecolor='#98c1a3', edgecolor='none', alpha=0.6)
    # water.plot(ax=ax1, facecolor='#a1cbe6', edgecolor='none', alpha=0.6)
    polygons.plot(ax=ax1, facecolor='none', edgecolor='red')

    # Annotation du GMC pour chaque polygone dans le premier graphique
    for idx, poly in polygons.iterrows():
        if 'GMC' in poly:
            x, y = poly.geometry.centroid.x, poly.geometry.centroid.y
            ax1.text(x, y, str(poly['GMC']), color='black', fontsize=8, ha='center')

    # # Personnalisation de l'arrière-plan
    ax1.set_facecolor('white')
    plt.show()

    # Graphique 2 : Rues, cours d'eau, espaces verts, avec polygones
    #fig2, ax2 = plt.subplots(figsize=(10, 10), facecolor='black',  dpi=100)
    #ax2.set_facecolor('black')
    #ox.plot_graph(streets, node_size=0, edge_color='white', bgcolor='black', edge_linewidth=0.8, show=False, close=False, ax=ax2)
    #green.plot(ax=ax2, facecolor='#98c1a3', edgecolor='none', alpha=0.6)
    #water.plot(ax=ax2, facecolor='#a1cbe6', edgecolor='none', alpha=0.6)
    #polygons.plot(ax=ax2, facecolor='none', edgecolor='white')

    # Annotation du GMC pour chaque polygone dans le second graphique
    # for idx, poly in polygons.iterrows():
    #     if 'GMC' in poly:
    #         x, y = poly.geometry.centroid.x, poly.geometry.centroid.y
    #         ax2.text(x, y, str(poly['GMC']), color='white', fontsize=8, ha='center')

    # Afficher les graphiques
    # ax2.set_facecolor('#a1cbe6')  # Bleu clair
    #plt.show()


#%% Copenhagen

coord_CPH = (55.6761, 12.5683)
CPH_path = r'C:\Users\Garance Dubuisson\OneDrive - Danmarks Tekniske Universitet\Thesis\QGIS\v4\CPH_districts2.shp'

stats_CPH = calculate_urban_parameters_from_shapefile(CPH_path)
# environ 45 min pour 60 polygones à CPH
stats_df = pd.DataFrame(stats_CPH) 

excel_path = r"C:\Users\Garance Dubuisson\OneDrive - Danmarks Tekniske Universitet\Thesis\CPH urban stat.xlsx"

#stats_df.to_excel(excel_path, sheet_name='v4', index=False)


with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a') as writer:
    workbook = writer.book
    if 'v4' in workbook.sheetnames:
        del workbook['v4']
    stats_df.to_excel(writer, sheet_name='v4', index=False)
    print("Les données ont été écrites dans la feuille 'v4'.")


get_pattern_polygons(CPH_path, coord_CPH)

#%% Work in progress

polygons = gpd.read_file(CPH_path).to_crs(epsg=4326)
#polygons = gpd.read_file(CPH_path).to_crs(epsg=3857)
#polygons = polygons.to_crs(epsg=25832)  # Utilisez le CRS UTM approprié pour votre région


for idx, poly in polygons.iterrows(): 
    polygon = poly.geometry
    #if polygon.is_valid:
        #print("Polygon is valid.")
      
    gmc = poly.get('GMC', None)
    print(f"Processing district with GMC: {gmc}")
    
    # in CRS 3857
    buildings = get_features_with_exception_from_polygon(polygon, tags_buildings)
    water = get_features_with_exception_from_polygon(polygon, tags_water)
    green = get_features_with_exception_from_polygon(polygon, tags_green)
    public = get_features_with_exception_from_polygon(polygon, tags_public)
    office = get_features_with_exception_from_polygon(polygon, tags_office)
    commercial = get_features_with_exception_from_polygon(polygon, tags_commercial)
            
    
    # Calculer les paramètres
    total_area = poly.get('Area_km2')*1000000 #in m2
    tot_area = polygon.area*1000000
    pop_dens = poly.get('Pop_densit')
    building_area = buildings.geometry.area.sum()
    GSI = building_area / total_area
    print("GSI from QGIS:", GSI, "vs with Python:", building_area/tot_area)


#%% Chat GPT
polygons = gpd.read_file(CPH_path).to_crs(epsg=3857)
polygons["tot_area"] = polygons.geometry.area

polygons = polygons.to_crs(epsg=4326)

for idx, poly in polygons.iterrows(): 
    polygon = poly.geometry
    #if polygon.is_valid:
        #print("Polygon is valid.")
      
    gmc = poly.get('GMC', None)
    print(f"Processing district with GMC: {gmc}")
    
    L_den = poly.get('Intersecti')
    L_night = poly.get('Intersec_1')
    AADT = poly.get('CPH_traf_1')
    if AADT is not None:
        AADT=AADT/total_area # veh/day/m²          
    Total_AADT = poly.get('CPH_traffi')
    print('noise levels:',L_den, L_night) 
    print('AADT:', AADT, Total_AADT)
