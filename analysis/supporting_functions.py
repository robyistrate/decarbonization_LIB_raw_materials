"""
Functions to create and manipulate LCI inventories and to perform LCIA calculations
"""

import numpy as np
import pandas as pd
import wurst
from constructive_geometries import *
import yaml
import copy
import brightway2 as bw
import presamples as ps
import geopandas as gpd
import pycountry


def get_db_as_dict(source_db: str):
    """
    Import the ecoinvent/biosphere database into a dictionary format
    """
    return [ds.as_dict() for ds in bw.Database(source_db)]


def get_breakdown_lists():
    with open("breakdown_lists.yaml", "r") as stream:
        try:
            data = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return data


def correct_product_in_exchanges(db):
    '''
    The function ExcelImporter requires that a 'reference product' is defined for every technosphere exchange.
    However, the ecoinvent database imported in BW2 uses 'product' for technosphere exchanges.
    
    This function modifies the datasets in place to replace 'reference product' by 'product' for
    all technosphere exchanges
    
    Arguments:
        - db: list of dictionaries; each dictionary is a dataset/activity
    '''    
    for ds in db:
        for exc in ds["exchanges"]:
            if 'reference product' in exc:
                exc['product'] = exc.pop('reference product')


def replicate_activity_to_location(ds, loc, db, DB_NAME):
    """
    Replicate an activity to a new location and relink exchanges to suppliers within the new location
    
    Arguments:
        ds (dict): The dataset that is replicated
        loc (str): The new location to replicate the activity to
        db (list): List of dicts containing the database activities for linking
        DB_NAME (str): Name of the database where the replicated dataset will be stored
    Returns:
        dictionary with the activity replicated to the new location.
    """
    production = lambda x: x["type"] == "production"

    # Replicate activity to the new locations
    ds_loc = wurst.transformations.geo.copy_to_new_location(ds, loc)

    # Translate the copy to the corresponding database
    ds_loc['database'] = DB_NAME

    # Change input code for production type
    for exc in filter(production, ds_loc["exchanges"]):
        exc.update({"database": ds_loc['database'],
                    'input': (ds_loc['database'], ds_loc['code'])})  
    
    # Relink exchanges:
    ds_loc_relink = relink_exchange_location(ds_loc, db)
 
    return ds_loc_relink


def get_dataset_for_location(exc_filter, loc, db):
    """
    Find dataset for a specific location

    :param exc_filter: dictionary containig the name, reference product, and unit of the searched dataset
    :param loc: location
    :param db: list of dictionaries containing the datasets
    """
    geomatcher = Geomatcher() # Initialize the geomatcher object    

    # Get the list of possible datasets for the exchange
    possible_datasets = list(wurst.transformations.geo.get_possibles(exc_filter, db)) 
        
    # Get both "market group" and "market" activities
    if 'market' in exc_filter['name']:
        exc_filter_market = copy.deepcopy(exc_filter)

        if "market for" in exc_filter['name']:
            exc_filter_market.update({'name': exc_filter['name'].replace('market', 'market group')})

        elif "market group" in exc_filter['name']:
            exc_filter_market.update({'name': exc_filter['name'].replace('market group', 'market')})

        possible_datasets = possible_datasets + list(wurst.transformations.geo.get_possibles(exc_filter_market, db))
        
    # Check if there is an exact match for the location
    match_dataset = [ds for ds in possible_datasets if ds['location'] == loc]

    if len(match_dataset) == 0:
        # If there is no specific dataset for the location, search for the supraregional locations
        loc_intersection = geomatcher.intersects(loc, biggest_first=False)
            
        for loc in [i[1] if type(i)==tuple else i for i in loc_intersection]:
            match_dataset = [ds for ds in possible_datasets if ds['location'] == loc]
            if len(match_dataset) > 0:
                break
            else:
                match_dataset = [ds for ds in possible_datasets if ds['location'] == 'RoW']
    return match_dataset[0]


def relink_exchange_location(ds, db):
    """
    Find new technosphere suppliers based on the location of the dataset.
    The new supplier is linked by code.

    Based on 'wurst.transformations.geo.relink_technosphere_exchanges'

    Arguments:
        ds (dict): The dataset
        dbs (List): List of datasets that contains the data to be used in the function
        
    Returns:
        The activity replicated to the new location.  
    """

    LOCATION = ds['location']
    technosphere = lambda x: x["type"] == "technosphere"
    geomatcher = Geomatcher() # Initialize the geomatcher object     

    for exc in filter(technosphere, ds["exchanges"]):
        exc_filter = {'name': exc['name'],
                      'product': exc['product'],
                      'unit': exc['unit']}

        match_dataset = get_dataset_for_location(exc_filter, LOCATION, db)

        exc.update({
                    'name': match_dataset['name'],
                    'product': match_dataset['reference product'],
                    'unit': match_dataset['unit'],
                    'location': match_dataset['location'],
                    'input': (match_dataset['database'], match_dataset['code'])
                    })
        
    # Include only unique exchanges (this is needed when the replicated dataset is RoW since it includes
    # many suppliers)

    unique_exch = []
    for exc in ds["exchanges"]:

        if exc["type"] == "technosphere":
            unique_exchanges_inputs = [i["input"] for i in unique_exch]
            if exc["input"] not in unique_exchanges_inputs:
                unique_exch.append(exc)
            else:
                exc_for_update = [e for e in unique_exch if e["input"] == exc["input"]][0]
                amount_update = exc_for_update["amount"] + exc["amount"]
                exc_for_update.update({'amount':amount_update,
                        })
        else:
            unique_exch.append(exc)

    ds["exchanges"] = unique_exch 
    return ds


def multi_lcia(activity, lcia_methods, amount=1):
    """
    Calculate multiple impact categories.
    
    Parameters:
    - activity (object): An activity object representing the product or process being assessed.
    - lcia_methods (dict): A dictionary of impact categories and their corresponding method. 
                           The keys are the names of the impact categories and the values are the methods to be used for each category.
    - amount (float): The functional unit of the assessment. Defaults to 1.
    
    Returns:
    - multi_lcia_results (dict): A dictionary of impact categories and their corresponding scores.
                                 The keys are the names of the impact categories and the values are the scores.
    """
    lca = bw.LCA({activity.key: amount})
    lca.lci()
    results = dict()
    for impact in lcia_methods:
        lca.switch_method(lcia_methods[impact])
        lca.lcia()
        results[impact] = lca.score
    return results


def lcia_direct_emissions(activity, lcia_methods, amount=1):
    method_CFs = {
        impact: {ef[0]: ef[1] 
                 for ef in bw.Method(lcia_methods[impact]).load()} 
                 for impact in lcia_methods
                 }

    emissions_impact = {
        impact: 0 for impact in lcia_methods
        }

    for exc in activity.biosphere():
        exc_amount = amount * exc['amount']
        exc_key = bw.get_activity(exc.input).key
        for impact in lcia_methods:
            if exc_key in method_CFs[impact]:
                emissions_impact[impact] += exc_amount * method_CFs[impact][exc_key]
    return emissions_impact


def lcia_heating_and_fuel(activity, lcia_methods, amount=1):

    fuels_list = [
        "natural gas, high pressure",
        "diesel",
        "hard coal",
        "light fuel oil"
    ]

    impacts_heating_and_fuel = {
        impact: {"Direct emissions": 0,
                 "Electricity": 0,
                 "Fuel": 0,
                 "Other": 0}
                 for impact in lcia_methods
                }

    # Direct emissons from assessed activity
    multi_lcia_results = lcia_direct_emissions(activity, lcia_methods, amount)
    for impact in multi_lcia_results:
        impacts_heating_and_fuel[impact]["Direct emissions"] += multi_lcia_results[impact]

    # ***********************************************
    if activity["reference product"] == "heat, from steam, in chemical industry":

        for exc in activity.technosphere():
            exc_amount = amount * exc['amount']

            if "heat, district or industrial" in exc.input["reference product"]:
                  
                # Direct emissions from each heat activity
                multi_lcia_results = lcia_direct_emissions(exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    impacts_heating_and_fuel[impact]["Direct emissions"] += multi_lcia_results[impact]

                # Emissions from inputs for each heat activity
                for eee in exc.input.technosphere():
                    eee_amount = exc_amount * eee['amount']
                    multi_lcia_results = multi_lcia(eee.input, lcia_methods, eee_amount)

                    # Electricity
                    if eee.input['reference product'] in ['electricity, low voltage', "electricity, medium voltage", "electricity, high voltage"]:
                        for impact in multi_lcia_results:
                            impacts_heating_and_fuel[impact]["Electricity"] += multi_lcia_results[impact]

                    # Fuel
                    elif eee.input['reference product'] in fuels_list:
                        for impact in multi_lcia_results:
                            impacts_heating_and_fuel[impact]["Fuel"] += multi_lcia_results[impact]

                    # Other
                    else:
                        for impact in multi_lcia_results:
                            impacts_heating_and_fuel[impact]["Other"] += multi_lcia_results[impact]

            # Emissions from non heat activities
            else:
                multi_lcia_results = multi_lcia(exc.input, lcia_methods, exc_amount)

                # Electricity
                if exc.input['reference product'] in ['electricity, low voltage', "electricity, medium voltage", "electricity, high voltage"]:
                    for impact in multi_lcia_results:
                        impacts_heating_and_fuel[impact]["Electricity"] += multi_lcia_results[impact]

                # Fuel
                elif exc.input['reference product'] in fuels_list:
                    for impact in multi_lcia_results:
                        impacts_heating_and_fuel[impact]["Fuel"] += multi_lcia_results[impact]

                # Other
                else:
                    for impact in multi_lcia_results:
                        impacts_heating_and_fuel[impact]["Other"] += multi_lcia_results[impact]

    # ***********************************************
    else:
    
        # Impact from inputs:
        for exc in activity.technosphere():
            exc_amount = amount * exc['amount']
            multi_lcia_results = multi_lcia(exc.input, lcia_methods, exc_amount)

            # Electricity
            if exc.input['reference product'] in ['electricity, low voltage',
                                                "electricity, medium voltage",
                                                "electricity, high voltage"]:
                for impact in multi_lcia_results:
                    impacts_heating_and_fuel[impact]["Electricity"] += multi_lcia_results[impact]

            # Fuel
            elif exc.input['reference product'] in fuels_list:
                for impact in multi_lcia_results:
                    impacts_heating_and_fuel[impact]["Fuel"] += multi_lcia_results[impact]

            # Other
            else:
                for impact in multi_lcia_results:
                    impacts_heating_and_fuel[impact]["Other"] += multi_lcia_results[impact]

    return impacts_heating_and_fuel


def lcia_system_contribution(activity, skip_inventories, lcia_methods, CONTRIBUTORS_LIST, breakdown_lists, activity_amount=1):
    '''
    This function computes the contribution of each system component to the total impact
    Based on: https://github.com/brightway-lca/brightway2/blob/master/notebooks/Contribution%20analysis%20and%20comparison.ipynb

    Parameters:
    - activity (object): An activity object representing the product or process being assessed.
    - lcia_methods (dict): A dictionary of impact categories and their corresponding method. 
                           The keys are the names of the impact categories and the values are the methods to be used for each category.
    - amount (float): The functional unit of the assessment. Defaults to 1.
    
    Returns:
    - system_contributions (dict): A nested dictionary of impact categories and system components and their corresponding LCIA scores.
                                   The keys are the names of the impact categories and the name of the system components, while the values are the scores.
    '''
   
    # Create am empty dict with the structure
    system_contributions = {
        impact: {
            contributor: 0 
            for contributor in CONTRIBUTORS_LIST} 
            for impact in lcia_methods
        }

    # Contribution of each system component
   
    # Process emissions
    process_emissions_impact = lcia_direct_emissions(activity, lcia_methods, amount=1)    
    for impact in process_emissions_impact:
        system_contributions[impact]['Direct emissions, process'] += process_emissions_impact[impact]

    for exc in activity.technosphere():
        
        # Skip inventories with breakdown to avoid double-counting
        if exc["name"] in skip_inventories:
            pass
        else:
            exc_amount = activity_amount * exc['amount']

            # Heating emissions
            if exc.input['reference product'] in breakdown_lists["heating products"]:
                heating_impacts = lcia_heating_and_fuel(exc.input, lcia_methods, exc_amount)
                for impact in heating_impacts:
                    system_contributions[impact]['Direct emissions, heat'] += heating_impacts[impact]["Direct emissions"]
                    system_contributions[impact]['Electricity consumption'] += heating_impacts[impact]["Electricity"]
                    system_contributions[impact]['Fuels consumption'] += heating_impacts[impact]["Fuel"]
                    system_contributions[impact]['Other'] += heating_impacts[impact]["Other"]

            # Fuels emissions
            elif exc.input['reference product'] in breakdown_lists["fuel products"]:
                fuel_impacts = lcia_heating_and_fuel(exc.input, lcia_methods, exc_amount)
                for impact in fuel_impacts:
                    system_contributions[impact]['Direct emissions, fuels'] += fuel_impacts[impact]["Direct emissions"]
                    system_contributions[impact]['Electricity consumption'] += fuel_impacts[impact]["Electricity"]
                    system_contributions[impact]['Fuels consumption'] += fuel_impacts[impact]["Fuel"]
                    system_contributions[impact]['Other'] += fuel_impacts[impact]["Other"]

            # Electricity
            elif exc.input['reference product'] in ['electricity, low voltage',
                                                "electricity, medium voltage",
                                                "electricity, high voltage",
                                                "diesel, burned in diesel-electric generating set, 10MW"]:
                multi_lcia_results = multi_lcia(exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    system_contributions[impact]['Electricity consumption'] += multi_lcia_results[impact]
            
            # Reagents
            elif exc.input['reference product'] in breakdown_lists["reagent products"]:        
                multi_lcia_results = multi_lcia(exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    system_contributions[impact]['Reagents consumption'] += multi_lcia_results[impact]      

            # Other activities
            else:       
                multi_lcia_results = multi_lcia(exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    system_contributions[impact]['Other'] += multi_lcia_results[impact]            

    return system_contributions


def lcia_reagents_disaggregation(activity, skip_inventories, lcia_methods, breakdown_lists, activity_amount=1):
    '''
    '''
   
    # Create am empty dict with the structure
    reagents_contributions = {
        impact: {
            reagent: 0 
            for reagent in breakdown_lists["reagent products"]} 
            for impact in lcia_methods
        }

    for exc in activity.technosphere():
        
        # Skip inventories with breakdown to avoid double-counting
        if exc["name"] in skip_inventories:
            pass
        else:
            exc_amount = activity_amount * exc['amount']
           
            if exc.input['reference product'] in breakdown_lists["reagent products"]:        
                multi_lcia_results = multi_lcia(exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    reagents_contributions[impact][exc.input['reference product']] = multi_lcia_results[impact]      
     

    return reagents_contributions