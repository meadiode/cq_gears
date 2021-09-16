
import pytest


def pytest_addoption(parser):
    parser.addoption('--rng_seed', action='store', type=int, default=42)
    parser.addoption('--num_test_cases', action='store', type=int, default=10)


def pytest_generate_tests(metafunc):
    funcarglist = metafunc.cls.gen_params(metafunc.config.option.rng_seed,
                                          metafunc.config.option.num_test_cases)

    argnames = metafunc.cls.argnames
    
    metafunc.parametrize(argnames,
        [[funcargs[name] for name in argnames] for funcargs in funcarglist])
