# -------------------------------------------------------------------------------
# Name: RKM Invasives Population Counter
# Purpose: Clusters species occurrences within 100m as "populations" and counts number of populations statewide and within PRISMs
# Author: Molly Moore
# Created: 2023-10-09
# Updated:
#
# To Do List/Future ideas:
# make more modular/universal in future if needed
# -------------------------------------------------------------------------------

# import packages and set environment variables
import arcpy
import os
import numpy as np
import datetime
import csv

arcpy.env.overwriteOutput = True
arcpy.env.workspace = "memory"

# set paths to data and working gdb - this stuff might change next time depending on input layers used and names of input layers
working_gdb = r"H:\\Projects\\RKM_Invasives\\Tiers\\Data\\final_layers.gdb"
prism = os.path.join(working_gdb, "PRISM_pa")
species_list = os.path.join(working_gdb, "PA_ultimate_species_list_2023_ids")

imap_points = os.path.join(working_gdb, r"iMap_Presence_Points_20230920")
imap_lines = os.path.join(working_gdb, r"iMap_Presence_Lines_20230920")
imap_polys = os.path.join(working_gdb, r"iMap_Presence_Polygons_20230920")
imap_approx = os.path.join(working_gdb, r"iMap_Approximate_Points_20230920")
gbif = os.path.join(working_gdb, r"gbif_export")
inat = os.path.join(working_gdb, r"inat_export")
nas = os.path.join(working_gdb, r"NAS_data_20230926_final")

lyr_list = [imap_points, imap_lines, imap_polys, imap_approx, gbif, inat, nas]
id_field = "iMap_ID_int"

def CalcPercentile(inputFeatureClass, percentileField, updateField):
    tbl_lyr = arcpy.MakeTableView_management(inputFeatureClass,"fc_lyr","{} IS NOT NULL".format(id_field))
    c_arr = arcpy.da.FeatureClassToNumPyArray(tbl_lyr, percentileField)
    arr = [float(x[0]) for x in np.ndarray.flatten(c_arr)]
    arr = [i for i in arr if i > 0]

    ## to create 3 rank for example
    p1 = np.percentile(arr, 33.33333, interpolation='lower')  # rank = 0
    p2 = np.percentile(arr, 66.66666, interpolation='lower')  # rank = 1
    p3 = np.percentile(arr, 100, interpolation='lower')  # rank = 2

    # use cursor to update the new rank field
    with arcpy.da.UpdateCursor(inputFeatureClass, [percentileField, updateField, id_field]) as cursor:
        for row in cursor:
            if row[2] is not None:
                if row[0] < p1 and row[0] > 0:
                    row[1] = 2  # rank 0
                elif p1 <= row[0] and row[0] < p2:
                    row[1] = 3
                elif p2 <= row[0] and row[0] <= p3:
                    row[1] = 4
                else:
                    pass
                cursor.updateRow(row)

# create paths to store temporary output buffered layers
lyr_out = []
for x in list(range(0, 6)):
    lyr_out.append(os.path.join("memory", "lyr" + str(x)))

# buffer input layers by 50m (50m is half the separation distance we are using) and dissolve if same species ID
for lyr, out in zip(lyr_list, lyr_out):
    print("Buffering... " + lyr)
    arcpy.Buffer_analysis(lyr, out, "50 Meters", dissolve_option="LIST", dissolve_field=id_field)

# set field map to drop unneeded fields - should probably just include iMap_ID_int because we can join other fields back at end
fieldMappings = arcpy.FieldMappings()

for lyr in lyr_out:
    fieldMappings.addTable(lyr)
fields = fieldMappings.fields
keep_fields = ['iMap_ID_int', 'PA_sciName', 'PA.comName', 'ITIS_TSN', 'iNat.Taxon_ID_NS', 'gbif_kyr']
for field in fieldMappings.fields:
    if field.name not in keep_fields:
        print("deleting... " + field.name)
        fieldMappings.removeFieldMap(fieldMappings.findFieldMapIndex(field.name))

# merge output buffered layers
data_merge = arcpy.Merge_management(lyr_out, os.path.join("memory", "data_merge"), fieldMappings)

# convert multipart polys to single part polygons
# population = arcpy.MultipartToSinglepart_management(data_merge, os.path.join(working_gdb, "population_clusters"))
population = os.path.join(working_gdb, "population_clusters")

# spatial join, summary statistics, and pivot table to find count of species population clusters in PRISMs
sp_join = arcpy.SpatialJoin_analysis(population, prism, os.path.join(working_gdb, "sp_join"), "JOIN_ONE_TO_MANY", "KEEP_COMMON", match_option="INTERSECT")
# sp_join = os.path.join(working_gdb,"sp_join")
sum_stat = arcpy.Statistics_analysis(sp_join, os.path.join(working_gdb, "summary_stats"), [[id_field, "Count"]], [id_field, "Region"])
pivot = arcpy.PivotTable_management(sum_stat, id_field, "Region", "FREQUENCY", os.path.join(working_gdb, "population_counts"))

# add field for total count in PA
arcpy.AddField_management(pivot, "PA_count", "LONG", "", "", "", "PA Count")

# get list of all regions
regions = sorted({row[0].replace(" ", "_") for row in arcpy.da.SearchCursor(prism, "Region")})
fieldCount = len(regions)
with arcpy.da.UpdateCursor(pivot, regions) as cursor:
    for row in cursor:
        for field in range(fieldCount):
            if row[field] is None:
                row[field] = 0
                cursor.updateRow(row)

# get total count of species population clusters in PA by adding all population numbers within PRISMs
fields = ["PA_count"] + regions
with arcpy.da.UpdateCursor(pivot, fields) as cursor:
    for row in cursor:
        row[0] = row[1] + row[2] + row[3] + row[4] + row[5] + row[6]
        cursor.updateRow(row)

for field in fields:
    arcpy.AddField_management(pivot, field + "_Tier", "SHORT")
    CalcPercentile(pivot, field, field+"_Tier")

# join other species information to pivot table
join_fields = ["PA_sciName", "PA_comName", "alt.name", "ITIS_TSN", "iNat_taxon_id_int", "gbif_taxon_key_int"]
arcpy.JoinField_management(pivot, id_field, species_list, id_field, join_fields)

arcpy.env.qualifiedFieldNames = False
species_lyr = arcpy.MakeTableView_management(species_list,"species_lyr")
pivot_lyr = arcpy.MakeTableView_management(pivot, "pivot_lyr")
arcpy.AddJoin_management(species_lyr,id_field,pivot_lyr,id_field)
species_lyr1 = arcpy.ExportTable_conversion(species_lyr,r"memory\\species_lyr1")

oid_fieldname = arcpy.Describe(species_lyr1).OIDFieldName
# create dictionary of all records
export_fields = [oid_fieldname, "iMap_ID_int", "PA_sciName", "PA_comName", "ITIS_TSN", "iNat_taxon_id_int", "gbif_taxon_key_int",
                 "PA_count", "PA_count_Tier", "North_Central", "North_Central_Tier", "North_East", "North_East_Tier", "North_West", "North_West_Tier", "South_Central", "South_Central_Tier", "South_East", "South_East_Tier", "South_West", "South_West_Tier"]
pivot_dict = {}
with arcpy.da.SearchCursor(species_lyr1,export_fields) as cursor:
    for row in cursor:
        pivot_dict[row[0]] = [row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9],row[10],row[11],row[12],row[13],row[14],row[15],row[16],row[17],row[18],row[19],row[20]]

# write to .csv file
with open(os.path.join(r"H:\Projects\RKM_Invasives\Tiers",'Tiers_'+time.strftime("%d%b%Y")+'.csv'), 'w', newline='') as csvfile:
    csv_output = csv.writer(csvfile)
    csv_output.writerow(export_fields)
    for key in sorted(pivot_dict.keys()):
        if key is not None:
            csv_output.writerow([key] + pivot_dict[key])
