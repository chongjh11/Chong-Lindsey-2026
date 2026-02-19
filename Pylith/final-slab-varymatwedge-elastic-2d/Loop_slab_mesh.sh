#!/bin/bash

# This script runs PyLith to get the different results due to different meshes
# Last modified on 15-Aug-2025 by Jeng Hann Chong
#
#   INSTRUCTIONS:
#   1. Change the CFG_TEMPLATE and MESH_LIST variables below as needed.
#   2. Change the material properties in the config file if needed within the "step.cfg" file.
#   3. Make sure the list of meshes have already been generated and stored in Textfiles/listmesh.txt
#   4. Run this script in a terminal: ./Loop_slab_mesh.sh (make sure it is executable: chmod +x Loop_slab_mesh.sh)

############################################################################

# User inputs (The config file and the list of meshes to be tested)
# Make sure the list of meshes have already been generated and stored in Textfiles/listmesh.txt
CFG_TEMPLATE="step_slab_greens_final.cfg"
MESH_LIST="Textfiles/listmesh.txt"

while IFS= read -r meshfile
do
    echo ""
    echo "-------------------------------------"    
    echo "Running PyLith with mesh: $meshfile .msh"
    echo "-------------------------------------" 
    echo ""
   
    # Extract basename without extension to use in temp config and output folder renaming
    meshname=$(basename "$meshfile" .msh)

    # Replace the line with the new mesh file
    sed -i '' "s/^reader.filename = .*\.msh$/reader.filename = .\/Mesh\/$meshfile.msh/" "$CFG_TEMPLATE" # > "$TMP_CFG"
    
    # Make a copy of the modified config file for clarity
    cat "$CFG_TEMPLATE" > "step_slab_greens_final_log.cfg"

    # Run PyLith
    pylith "step_slab_greens_final_log.cfg"

    # Change the output directory name & move the temp config file into the output directory
    mv elastic_slab "elastic_slab_${meshname}"
    mv step_slab_greens_final_log.cfg "elastic_slab_${meshname}"

done < "$MESH_LIST"
