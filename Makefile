all:
	echo "Building distribution"
	python setup.py build

develop:
	python setup.py develop

clean:
	@rm -f jip/*.so
	@rm -f jip/*.pyc
	@rm -f *.pyc
	@rm -rf build/
