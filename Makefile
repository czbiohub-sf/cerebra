help:
	@echo "lint - check code style with flake8"
	@echo "test - run tests only"
	@echo "coverage - run tests and check code coverage"
	@echo "conda_install (recommended) - Install requirements "

test:
	py.test --ignore=test_mutation_counts.py

coverage:
	coverage run --source cerebra --omit="*/test*" --module py.test
	coverage report --show-missing

lint:
	flake8 --exclude docs cerebra --exit-zero # need to get rid of this

conda_install:
	conda install --file conda_requirements.txt
	pip install -r requirements.txt
	pip install .

install:
	pip install -r requirements.txt
	pip install .
