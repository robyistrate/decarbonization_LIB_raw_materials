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
import math


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


def link_exchanges_by_code(db, external_db, biosphere_db):
    '''
    This function links in place technosphere exchanges within the database and/or to an external database
    and biosphere exchanges with the biosphere database (only unlinked exchanges)
    
    Returns a dictionary with the linked database
    '''   
    technosphere = lambda x: x["type"] == "technosphere"
    biosphere = lambda x: x["type"] == "biosphere"
    
    for ds in db:
        
        for exc in filter(technosphere, ds["exchanges"]):
            if 'input' not in exc:
                try:
                    exc_lci = wurst.get_one(db + external_db,
                                            wurst.equals("name", exc['name']),
                                            wurst.equals("reference product", exc['product']),
                                            wurst.equals("location", exc['location'])
                                        )
                    exc.update({'input': (exc_lci['database'],
                                        exc_lci['code'])})
                except Exception:
                    print(exc['name'], exc['product'], exc['location'])
                    raise
            
        for exc in filter(biosphere, ds["exchanges"]):
            if 'input' not in exc:
                try:
                    ef_code = [ef['code'] for ef in biosphere_db if ef['name'] == exc['name'] and 
                                                                    ef['unit'] == exc['unit'] and 
                                                                    ef['categories'] == exc['categories']][0]
                    exc.update({'input': ('biosphere3',
                                          ef_code)})   
                except Exception:
                    print(exc['name'], exc['unit'], exc['categories'])
                    raise


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
        
        # There are no regions intersecting with "RoW"; check for "RoW" or "GLO datasets"
        if len(loc_intersection) == 0:
            loc_intersection = ["RoW", "GLO"]
            
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


def create_renewable_electricity_country_market(renewable_tech, countries, dbs, db_name):
    """
    This function creates country markets for renewable power generation technologies
    for those countries that have additional regionalization
    """
    technosphere = lambda x: x["type"] == "technosphere"
    renewable_tech_country_markets = []

    for c in countries:

        # Create market activity for country
        actv_dict = {
            "name": "market for " + renewable_tech,
            "reference product": "electricity, high voltage",
            "location": c,
            "unit": "kilowatt hour",
            "database": db_name,
            "code": wurst.filesystem.get_uuid()
        }
        
        exchanges = []

        # Add production to exchanges
        exchanges.append({
            'name': actv_dict["name"],
            'product': actv_dict['reference product'],
            "location": actv_dict['location'],
            'amount': 1,
            'unit': actv_dict['unit'],
            'database': actv_dict['database'],
            'type': "production"
            })
                    
        # Get market for electricity high voltage in the country
        elec_market_country = wurst.get_one(
            dbs,
            wurst.equals("name", "market group for electricity, high voltage"),
            wurst.equals("reference product", "electricity, high voltage"),
            wurst.equals("location", c)
            )

        # Iterate over regional markets within the country
        # and get market share for the region and associated dataset
        for reg_elec in elec_market_country["exchanges"]:
            elec_market_reg = wurst.get_one(
                dbs,
                wurst.equals("name", reg_elec["name"]),
                wurst.equals("reference product", reg_elec["product"]),
                wurst.equals("location", reg_elec["location"])
                )
            
            # Some markets have a second layer of regionalized markets
            if reg_elec["name"] == "market group for electricity, high voltage":
                for reg_elec_second in elec_market_reg["exchanges"]:
                    elec_market_reg_second = wurst.get_one(
                        dbs,
                        wurst.equals("name", reg_elec_second["name"]),
                        wurst.equals("reference product", reg_elec_second["product"]),
                        wurst.equals("location", reg_elec_second["location"])
                        )
                    
                    # Iterate over power generation technologies within
                    # each reginoal market
                    for tech_elec in elec_market_reg_second["exchanges"]:
                        if tech_elec["name"] == renewable_tech:
                            elec_share_tech = tech_elec["amount"] * reg_elec_second["amount"] * reg_elec["amount"]

                            exchanges.append({'name': tech_elec["name"],
                                            'product': tech_elec['product'],
                                            'amount': elec_share_tech,
                                            'unit': tech_elec['unit'],
                                            'database': tech_elec['database'],
                                            'location': tech_elec['location'],
                                            'type': tech_elec['type'],
                                            })

            # Other markets do not have further layers   
            # Iterate over power generation technologies within
            # each reginoal market
            for tech_elec in elec_market_reg["exchanges"]:
                
                if tech_elec["name"] == renewable_tech:
                    elec_share_tech = tech_elec["amount"] * reg_elec["amount"]

                    exchanges.append({'name': tech_elec["name"],
                                      'product': tech_elec['product'],
                                      'amount': elec_share_tech,
                                      'unit': tech_elec['unit'],
                                      'database': tech_elec['database'],
                                      'location': tech_elec['location'],
                                      'type': tech_elec['type']
                                      }) 
                    
        actv_dict.update({'exchanges': exchanges})
        renewable_tech_country_markets.append(actv_dict)

    # Normalize market shares:
    for country in renewable_tech_country_markets:
        # Get sum of technology shares
        sum_tech_shares = 0
        for exc in filter(technosphere, country["exchanges"]):
            sum_tech_shares += exc["amount"]
        
        # Normalize:
        for exc in filter(technosphere, country["exchanges"]):
            exc.update({
                "amount": exc["amount"] / sum_tech_shares
            })

    # Link exchanges:
    link_exchanges_by_code(renewable_tech_country_markets, dbs, [])

    return renewable_tech_country_markets


def create_pedigree_matrix(pedigree_scores: tuple, exc_amount: float):
    """
    This function returns a dict containing the pedigree matrix dict and loc and scale values
    that can be used to update exchanges in a dataset dict
    
    The pedigree matrix dictionary is created using the scores provided in the LCI Excel file.

    The code to calcualte the loc and scale values is based on https://github.com/brightway-lca/pedigree_matrix,
    which is published by Chris Mutel under an BSD 3-Clause License (2021).

    
    :param pedigree_scores: tuple of pedigree scores
    :param exc_amount: exchange amount
    :return dict:
    """

    VERSION_2 = {
        "reliability": (1.0, 1.54, 1.61, 1.69, 1.69),
        "completeness": (1.0, 1.03, 1.04, 1.08, 1.08),
        "temporal correlation": (1.0, 1.03, 1.1, 1.19, 1.29),
        "geographical correlation": (1.0, 1.04, 1.08, 1.11, 1.11),
        "further technological correlation": (1.0, 1.18, 1.65, 2.08, 2.8),
        "sample size": (1.0, 1.0, 1.0, 1.0, 1.0),
    }
    
    pedigree_scores_dict = {
        'reliability': pedigree_scores[0],
        'completeness': pedigree_scores[1],
        'temporal correlation': pedigree_scores[2],
        'geographical correlation': pedigree_scores[3],
        'further technological correlation': pedigree_scores[4]
    }
    
    assert len(pedigree_scores) in (5, 6), "Must provide either 5 or 6 factors"
    if len(pedigree_scores) == 5:
        pedigree_scores = pedigree_scores + (1,)

    factors = [VERSION_2[key][index - 1] for key, index in pedigree_scores_dict.items()]

    basic_uncertainty: float = 1.0
    values = [basic_uncertainty] + factors

    scale = math.sqrt(sum([math.log(x) ** 2 for x in values])) / 2
    loc = math.log(abs(exc_amount))

    pedigree_dict = {
        'uncertainty type': 2,
        'loc': loc,
        'scale': scale,
        "pedigree": pedigree_scores_dict,
    }
    return pedigree_dict


def init_simple_lca(activity):
    """
    Initialize simple LCA object
    """
    lca = bw.LCA({activity.key: 1})
    lca.lci()
    return lca


def init_mc_lca(activity):
    """
    Initialize Monte Carlo LCA object

    Create a MonteCarloLCA object with functional unit but no method. 
    Run .lci to build the A and B matrices (and hence fix the indices of our matrices)
    """
    mc_lca = bw.MonteCarloLCA({activity: 1})
    mc_lca.lci()
    return mc_lca


def multi_lcia(lca, activity, lcia_methods, amount=1):
    """
    Calculate multiple impact categories.
    
    Parameters:
    - lca: lca object
    - activity (object): An activity object representing the product or process being assessed.
    - lcia_methods (dict): A dictionary of impact categories and their corresponding method. 
                           The keys are the names of the impact categories and the values are the methods to be used for each category.
    - amount (float): The functional unit of the assessment. Defaults to 1.
    
    Returns:
    - multi_lcia_results (dict): A dictionary of impact categories and their corresponding scores.
                                 The keys are the names of the impact categories and the values are the scores.
    """

    results = dict()
    lca.redo_lci({activity.key: amount})

    for impact in lcia_methods:
        lca.switch_method(lcia_methods[impact])
        lca.lcia()
        results[impact] = lca.score
    return results

def multi_lcia_MonteCarlo(mc_lca, activity, iterations, lcia_methods, amount=1):
    '''
    This function performs Monte Carlo simulation accross a range of categories

    :param activity: dataset of the activity for analysis
    :param iterations: number of Monte Carlo runs
    :param lcia_methods: dictionary containing assessed impact categories and methods
    '''

    mc_lca.redo_lci({activity.key: amount})
    
    # Iterate over the methods and stores the C matrix in a dictionary
    C_matrices = dict()
    for method in lcia_methods:
        mc_lca.switch_method(lcia_methods[method])
        C_matrices[method] = mc_lca.characterization_matrix
        
    # Create a dictionary to store the results
    results_MC = np.empty((len(lcia_methods), iterations))
   
    # Populate results dictionary using a for loop over the number of MonteCarlo iterations required 
    # Include a nested for loop over the methods and multiply the characterization matrix with the inventory.       
    for iteration in range(iterations):
        next(mc_lca)
        for method_index, method in enumerate(lcia_methods):
            results_MC[method_index, iteration] = (C_matrices[method] * mc_lca.inventory).sum()
       
    results_MC_dict = dict()
    i_count = 0
    for method in lcia_methods:
        results_MC_dict[method] = results_MC[i_count]
        i_count += 1

    return results_MC_dict


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


def lcia_system_contribution(lca, activity, skip_inventories, lcia_methods, CONTRIBUTORS_LIST, breakdown_lists, activity_amount=1):
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
    process_emissions_impact = lcia_direct_emissions(activity, lcia_methods, activity_amount)    
    for impact in process_emissions_impact:

        # Nickel datasets require some special treatment:

        # For the nickel mining dataset,
        # direct emissions are attributed to "Fuels consumption":
        if activity["name"] == "nickel ore mining, average excluding China":
            system_contributions[impact]['Fuels consumption'] += process_emissions_impact[impact]
        
        # For nickel concentration dataset,
        # 15% of direct emissions are attributed to "Fuels consumption"
        # and 85% to "Process emissions" as these are emissions from the use of flotation agents containing C
        elif activity["name"] == "nickel concentration, average excluding China":
            system_contributions[impact]['Fuels consumption'] += process_emissions_impact[impact] * 0.15
            system_contributions[impact]['Process emissions'] += process_emissions_impact[impact] * 0.85
    
        # For nickel matte production,
        # 3% of direct emissions are attributed to "Fuels consumption"
        # and 97% to "Process emissions"
        elif activity["name"] == "nickel matte production, nickel-cobalt sulphite and nickel sub materials, mass allocation":
            system_contributions[impact]['Fuels consumption'] += process_emissions_impact[impact] * 0.03
            system_contributions[impact]['Process emissions'] += process_emissions_impact[impact] * 0.97            
        
        # For nickel sulfate production, all direct emissions are process emissions

        else:
            system_contributions[impact]['Process emissions'] += process_emissions_impact[impact]

    for exc in activity.technosphere():
        
        # Skip inventories with breakdown to avoid double-counting
        if exc["name"] in skip_inventories:
            pass
        else:
            exc_amount = activity_amount * exc['amount']

            # Electricity consumption
            if exc.input['reference product'] in breakdown_lists["electricity products"]:
                multi_lcia_results = multi_lcia(lca, exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    system_contributions[impact]['Electricity consumption'] += multi_lcia_results[impact]
                
            # Process heating emissions
            elif exc.input['reference product'] in breakdown_lists["heating products"]:
                multi_lcia_results = multi_lcia(lca, exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    system_contributions[impact]['Process heating'] += multi_lcia_results[impact]            

            # Fuels consumption
            elif exc.input['reference product'] in breakdown_lists["fuel products"]:
                multi_lcia_results = multi_lcia(lca, exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    system_contributions[impact]['Fuels consumption'] += multi_lcia_results[impact]  

            # Reagents consumption
            elif exc.input['reference product'] in breakdown_lists["reagent products"]:        
                multi_lcia_results = multi_lcia(lca, exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    system_contributions[impact]['Reagents consumption'] += multi_lcia_results[impact]    

            # Other activities
            else:       
                multi_lcia_results = multi_lcia(lca, exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    system_contributions[impact]['Other'] += multi_lcia_results[impact]   

    return system_contributions


def lcia_reagents_disaggregation(lca, activity, skip_inventories, lcia_methods, breakdown_lists, activity_amount=1):
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
                multi_lcia_results = multi_lcia(lca, exc.input, lcia_methods, exc_amount)
                for impact in multi_lcia_results:
                    reagents_contributions[impact][exc.input['reference product']] += multi_lcia_results[impact]      
     

    return reagents_contributions