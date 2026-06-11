import pybdsim
import os


def test():

    os.chdir(os.path.dirname(__file__))

    # === Setup filenames ===
    base_name = "field_1d"
    template_name = "field_track.tpl"
    gmad_name = base_name + ".gmad"
    root_name = base_name + ".root"

    # === Test parameters ===
    data = {
        'FIELD_TYPE': 'bmap',
        'FIELD_FILE': '1dexample.dat',
        'FIELD_LENGTH': '1.0',
        'BEAM_ENERGY': '1.0'
    }

    pybdsim.Run.RenderGmadJinjaTemplate(template_name, gmad_name, data)
    pybdsim.Run.Bdsim(gmad_name, base_name, ngenerate=5000, seed=1)

    d = pybdsim.DataPandas.BDSIMOutput(root_name)