#!/bin/bash

# This script runs PyLith to get the different results due to different meshes
# Last modified on 5-Nov-2025 by Jeng Hann Chong
#
#   INSTRUCTIONS:
#   1. Change the CFG_TEMPLATE and SPD_LIST variables below as needed.
#   2. Change the material properties in the config file if needed within the "step.cfg" file.
#   3. Make sure the list of meshes have already been generated and stored in Textfiles/listmatwed.txt
#   4. Run this script in a terminal: ./Loop_slab_sig0_mesh.sh (make sure it is executable: chmod +x Loop_slab_mesh.sh)

############################################################################

# User inputs (The config file and the list of meshes to be tested)
# Make sure the list of meshes have already been generated and stored in Textfiles/listmesh.txt
#!/bin/bash

CFG_TEMPLATE="step_slab_greens_sig0_final.cfg"
SPD_LIST="Textfiles/listmatwed20deg.txt"

while IFS= read -r spatialdbfile; do
    spatialdbname=$(basename "$spatialdbfile" .spatialdb)
    TMP_CFG="temp_${spatialdbname}.cfg"
    cp "$CFG_TEMPLATE" "$TMP_CFG"

    echo ""
    echo "-------------------------------------"    
    echo "Running PyLith with varying prism property: $spatialdbfile"
    echo "-------------------------------------" 
    echo ""

    # Replace only the line inside continental_mantle section
    # CHANGE ACCORDING to the path if needed at the "awk" line below
    awk -v file="./ShearMods/20deg40thc/$spatialdbfile" '

    BEGIN { in_section=0 }
    /^\[pylithapp\.problem\.materials\.continental_mantle\]/ { in_section=1 }
    /^\[pylithapp\.problem\.materials\./ && !/\[pylithapp\.problem\.materials\.continental_mantle\]/ { in_section=0 }
    in_section && /^db_auxiliary_field\.iohandler\.filename =/ {
        print "db_auxiliary_field.iohandler.filename = " file
        next
    }
    { print }
    ' "$TMP_CFG" > "${TMP_CFG}.tmp" && mv "${TMP_CFG}.tmp" "$TMP_CFG"

    # Run PyLith
    pylith "$TMP_CFG"

    # Rename output folder and move the modified config inside
    if [ -d elastic_slab_sig0 ]; then
        mv elastic_slab_sig0 "elastic_slab_sig0_${spatialdbname}"
    fi
    mv "$TMP_CFG" "./elastic_slab_sig0_${spatialdbname}/"

    echo ""
    echo "-------------------------------------"    
    echo "Completed: $spatialdbfile"
    echo "-------------------------------------" 
    echo ""

done < "$SPD_LIST"
