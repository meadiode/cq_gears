import pytest
import json


def pytest_generate_tests(metafunc):

    cases = json.loads(open('regression_test_cases.json').read())

    metafunc.parametrize(('cls_name', 'expected', 'args',),
                         [(case['class'],
                           case['expected'],
                           case['args'],) for case in cases])
