install:
	python setup.py install

dev:
	python setup.py develop

deps:
	pip install -r requirements.txt

inplace:
	python setup.py build_ext --inplace

test:
	pytest tests

clean:
	python setup.py clean --all
	$(RM) ./distruct/_diSTruct.cpp _diSTruct.*.so
	rm -rf ./shallonwk

.PHONY: install dev deps inplace test clean
