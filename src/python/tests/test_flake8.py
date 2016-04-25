import subprocess
import os

def test_flake8():
    python_grackle_dir = os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))
    output_file = os.environ.get("WORKSPACE", None) or os.getcwd()
    output_file = os.path.join(output_file, 'flake8.out')
    if os.path.exists(output_file):
        os.remove(output_file)
    output_string = "--output-file=%s" % output_file
    config_string = "--config=%s" % os.path.join(
        os.path.dirname(python_grackle_dir), 'setup.cfg')
    subprocess.call(['flake8', output_string, config_string, python_grackle_dir])

    with open(output_file) as f:
        flake8_output = f.readlines()
    if flake8_output != []:
        raise AssertionError(
            "flake8 found style errors:\n\n%s" % "\n".join(flake8_output))
