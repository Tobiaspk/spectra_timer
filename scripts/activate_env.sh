# if venv does not exist, create it
if [ ! -d "venv" ]; then
    python3 -m venv venv
fi

# if venv is not activated, activate it
if [ -z "$VIRTUAL_ENV" ]; then
    source venv/bin/activate
fi

# install and upgrade this package
pip install --upgrade .

# install line_profiler
pip install line_profiler