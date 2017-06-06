set PREFIX = '/Users/Giuliano/anaconda/envs'

$PREFIX/py33/bin/./python setup.py bdist_wheel --universal
$PREFIX/py33/bin/./python setup.py sdist
twine upload dist/*
rm -r build
rm -r roteas.egg-info
rm -r dist

