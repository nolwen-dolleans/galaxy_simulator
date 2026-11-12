source venv/bin/activate

pip install -U pip
pip install -r requirements.txt

cmake --build build

./build/main $1

python script/galaxy_visualisation.py