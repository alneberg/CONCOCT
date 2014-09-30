.PHONY: clean-pyc clean-build docs clean

help:
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests with nosetests"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	@echo "release - package and upload a release"
	@echo "dist - package"

clean: clean-build clean-pyc
	@echo "clean!"

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

lint:
	flake8 concoct bin tests

test:
	nosetests

docs:
	python doc/add_cli_arguments_to_docs.py
	$(MAKE) -C doc clean
	$(MAKE) -C doc html
	open doc/build/html/index.html

release: clean
	python setup.py sdist upload

#docker:
#python doc/generate_dockerfile.py doc/Dockerfile.template doc/Dockerfile <version>
#boot2docker up
#scp doc/Dockerfile docker@192.168.59.103:~/Dockerfile
#boot2docker ssh
#docker build .
#docker images
#docker tag <docker id> binpro/concoct:<version>

dist: clean
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist
