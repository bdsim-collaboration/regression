import pytest
import os

pytestmark = pytest.mark.xfail(reason="requires bdsim")

def func_test():
    print("func_test")

def test_subprocess(make_bdsim_test_code, run_bdsim_test_code_as_subprocess):

    code_to_run = make_bdsim_test_code(func_test, args="", dir=os.path.dirname(os.path.abspath(__file__)))
    output = run_bdsim_test_code_as_subprocess(code_to_run)



