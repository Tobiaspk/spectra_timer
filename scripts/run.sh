# activate environment and reload package
zsh scripts/activate_env.sh
USE_CELL_TYPES=False

# run data preparations
python scripts/prep_data.py --use_cell_types $USE_CELL_TYPES 

# run the model
python scripts/run_model.py --use_cell_types $USE_CELL_TYPES  -e 1

# install pip in developer
if [ $USE_CELL_TYPES=="True" ];
then
    performance_path="scripts/performance_results_cell_types.txt"
else
    performance_path="scripts/performance_results.txt"
fi
kernprof -l -v scripts/run_model.py --use_cell_types $USE_CELL_TYPES  -e 100 > $performance_path
