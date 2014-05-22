all:
	echo "Building distribution"
	python setup.py build

develop:
	python setup.py develop

tests:
	py.test

mysqltest:
	py.test -m mysqltest --mysql "mysql:///test"

doc:
	$(MAKE) -C docs html

clean:
	@rm -f jip/*.so
	@rm -f jip/*.pyc
	@rm -f *.pyc
	@rm -rf build/
	@rm -rf dist/
	@rm -rf jip.egg-info/
	@rm -rf pyjip.egg-info/
	$(MAKE) -C docs clean
