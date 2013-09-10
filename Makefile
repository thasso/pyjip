all:
	echo "Building distribution"
	python setup.py build

develop:
	python setup.py develop

doc:
	$(MAKE) -C docs html

clean:
	@rm -f jip/*.so
	@rm -f jip/*.pyc
	@rm -f *.pyc
	@rm -rf build/
	@rm -rf jip.egg-info/
	$(MAKE) -C docs clean
