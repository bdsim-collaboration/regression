import os
from itertools import product

import numpy as np
import pandas as pd
import pybdsim
import pytest
from pybdsim.Testing.RotateMatrix import (
    compose_order2,
    edge_R_matrix,
    edge_T_matrix,
    roll_matrix_bdsim,
    rotate_order2_map_lab_to_magnet,
    sector_body_T_minimal,
    rotated_matrix
)

from pybdsim.Testing.CompareMatrix import (
    compare_matrix,
    round_matrix,
    max_matrix_diff,
)


REGRESSION_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = REGRESSION_DIR
SIMULATION_RESULTS_CSV = os.path.join(RESULTS_DIR, "sbend_simulation_results.csv")
ORDER2_QUAD_PAIRS = (
    (0, 0),
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 2),
    (2, 3),
    (3, 3),
)

mass = {"e-": 0.51099895, "proton": 938.272}  
charge = {"e-": -1, "proton": 1}  

INTEGRATOR_PRECISIONS = {
    "bdsimmatrix": { #default: BDSIM dipole matrix
        "linear_precision": 1.0e-4,
        "linear_precision_tilt": 1.0e-6,
        "linear_precision_body": 1.0e-4,
        "order2_precision": 1.0e-2,
        "linear_precision_fint": 1.0e-5,
        "order2_precision_fint": 0.13,
        "order2_precision_fint2": 1.0e-1,
        "order2_precision_h1": 1.14,
    },
    "bdsim": { #BDSIM Dipole Rodrigues 2
        "linear_precision": 1.0e-4,
        "linear_precision_tilt": 1.0e-6,
        "linear_precision_body": 1.0e-4,
        "order2_precision": 1.0e-2,
        "linear_precision_fint": 1.0e-5,
        "order2_precision_fint": 1.0e-2,
        "order2_precision_fint2": 1.0e-1,
        "order2_precision_h1": 0.185,
    },
    "geant4": { #G4ClassicalRK4
        "linear_precision": 0.15,
        "linear_precision_tilt": 1.0e-4,
        "linear_precision_body": 1.0e-4,
        "order2_precision": 1.0e-1,
        "linear_precision_fint": 1.0e-6,
        "order2_precision_fint": 1.0e-2,
        "order2_precision_fint2": 1.0e-1,
        "order2_precision_h1": 0.43,
    },
    "bdsimtwo": { #acceptable if no fringe dipole pole faces are present. Dipole: BDSIM Dipole Rodrigues 2 and Dipole with K1: BDSIM Dipole Matrix
        "linear_precision": 0.15,
        "linear_precision_tilt": 1.0e-4,
        "linear_precision_body": 1.0e-4,
        "order2_precision": 1.0e-1,
        "linear_precision_fint": 1.0e-6,
        "order2_precision_fint": 1.0e-2,
        "order2_precision_fint2": 1.0e-1,
        "order2_precision_h1": 0.43,
    },
    "bdsimmatrixfringescaling": { # BDSIM dipole matrix
        "linear_precision": 1.0e-4,
        "linear_precision_tilt": 1.0e-6,
        "linear_precision_body": 1.0e-4,
        "order2_precision": 1.0e-2,
        "linear_precision_fint": 1.0e-4,
        "order2_precision_fint": 0.13,
        "order2_precision_fint2": 1.0e-1,
        "order2_precision_h1": 1.14,
    },
    "geant4dp": { #G4DormandPrince745
        "linear_precision": 0.15,
        "linear_precision_tilt": 1.0e-4,
        "linear_precision_body": 1.0e-4,
        "order2_precision": 1.0e-1,
        "linear_precision_fint": 1.0e-6,
        "order2_precision_fint": 1.0e-2,
        "order2_precision_fint2": 1.0e-1,
        "order2_precision_h1": 0.43,
    },
}


def make_config_id(config):
    angle_deg = np.degrees(config["ANGLE"])
    return (
        f'{config["PARTICLE_TYPE"]}_'
        f'{config["INTEGRATOR"]}_'
        f'L{config["LENGTH"]:g}_'
        f'A{angle_deg:g}_'
        f'E{config["BEAM_ENERGY"]:g}'
    )

def build_sbend_params_dataframe(
    lengths,
    angles,
    beam_energies,
    particles,
    integrators,
    tilt1,
    e1,
    e2,
    h1,
    h2,
    fint,
    fintx,
    fintk2,
    fintxk2,
    hgap,
    return_dict=False,
):

    if abs(e1) > np.pi / 4.0 or abs(e2) > np.pi / 4.0:
        print("Warning: pole face angles greater than 45 degrees may cause geometry issues")

    rows = []
    for particle_type, integrator, length, angle, beam_energy in product(
        particles, integrators, lengths, angles, beam_energies
    ):
        rho = length / angle
        horizontal_width = min(0.6, rho / 2.0 - 1.0e-3)
        rows.append(
            {
                "LENGTH": length,
                "ANGLE": angle,
                "BEAM_ENERGY": beam_energy,
                "PARTICLE_TYPE": particle_type,
                "HORIZONTAL_WIDTH": horizontal_width,
                "TILT1": tilt1,
                "E1": e1,
                "E2": e2,
                "INTEGRATOR": integrator,
                "H1": h1,
                "H2": h2,
                "FINT": fint,
                "FINTX": fintx,
                "FINTK2": fintk2,
                "FINTXK2": fintxk2,
                "HGAP": hgap,
                "BODY_K1": 0.0,
                **INTEGRATOR_PRECISIONS[integrator],
            }
        )
    df = pd.DataFrame(rows)
    df["first_order_body_reference"] = df["linear_precision_body"]
    df["first_order_reference"] = df["linear_precision"]
    df["second_order_tilt_linear_reference"] = df["linear_precision_tilt"]
    df["second_order_tilt_order2_reference"] = df["order2_precision"]
    df["h1_second_order_reference"] = df["order2_precision_h1"]
    df["fint_second_order_reference"] = df["order2_precision_fint"]
    df["fint_k2_second_order_reference"] = df["order2_precision_fint2"]
    df["second_order_fint_linear_reference"] = df["linear_precision_fint"]
    df["config_row_index"] = np.arange(len(df))
    df["config_id"] = [make_config_id(row) for row in df.to_dict("records")]
    
    if return_dict:
        return df, df.to_dict("list")
    return df


def build_sbend_test_configs():
    return build_sbend_params_dataframe(
        lengths=np.array([ 0.5, 1.5]), #0.5
        angles=np.radians(np.array([10.0, 40.0])), #10 ,40,# degrees
        beam_energies=np.array([10.0]),
        particles=["e-", "proton"],
        integrators=[
            "bdsimmatrix",
            "bdsimtwo",
            "bdsimmatrixfringescaling",
            "geant4",
            "geant4dp",
        ],
        tilt1=np.pi / 6.0,
        e1=np.pi / 15.0,
        e2=-np.pi / 15.0,
        h1=0.5,
        h2=0,
        fint=0.5,
        fintx=0.0,
        fintk2=0.5,
        fintxk2=0.5,
        hgap=0.06,
    )


CONFIG_DF = build_sbend_test_configs()
CONFIGS = CONFIG_DF.to_dict("records")
simulation_results = []
SIMULATION_RESULTS_CSV = os.path.join(REGRESSION_DIR, "Results", "sbend_results.csv")


def append_simulation_result(config, test_name, measured, reference):
    row = dict(config)          
    row["test_name"] = test_name
    row["measured"] = measured
    row["reference"] = reference
    row["passed"] = measured <= reference
    simulation_results.append(row)


def write_simulation_results_csv():
    os.makedirs(os.path.dirname(SIMULATION_RESULTS_CSV), exist_ok=True)
    pd.DataFrame(simulation_results).to_csv(
        SIMULATION_RESULTS_CSV, index=False, float_format="%.4g"
    )


@pytest.fixture(autouse=True)
def save_simulation_results():
    yield
    write_simulation_results_csv()

def CalculateTaylorMapOrder2(root_file_name, sampler1_name, sampler2_name, average=True):
    root_file = pybdsim.Data.Load(root_file_name)
    

    sampler1_data = pybdsim.Data.SamplerData(root_file, sampler1_name)
    sampler2_data = pybdsim.Data.SamplerData(root_file, sampler2_name)

    sampler1_x = np.asarray(sampler1_data.data["x"], dtype=float)
    sampler1_xp = np.asarray(sampler1_data.data["xp"], dtype=float)
    sampler1_y = np.asarray(sampler1_data.data["y"], dtype=float)
    sampler1_yp = np.asarray(sampler1_data.data["yp"], dtype=float)

    sampler2_x = np.asarray(sampler2_data.data["x"], dtype=float)
    sampler2_xp = np.asarray(sampler2_data.data["xp"], dtype=float)
    sampler2_y = np.asarray(sampler2_data.data["y"], dtype=float)
    sampler2_yp = np.asarray(sampler2_data.data["yp"], dtype=float)

    if average:
        sampler1_x = sampler1_x - sampler1_x.mean()
        sampler1_xp = sampler1_xp - sampler1_xp.mean()
        sampler1_y = sampler1_y - sampler1_y.mean()
        sampler1_yp = sampler1_yp - sampler1_yp.mean()

        sampler2_x = sampler2_x - sampler2_x.mean()
        sampler2_xp = sampler2_xp - sampler2_xp.mean()
        sampler2_y = sampler2_y - sampler2_y.mean()
        sampler2_yp = sampler2_yp - sampler2_yp.mean()

    sampler1_matrix = np.array(
        [
            sampler1_x,
            sampler1_xp,
            sampler1_y,
            sampler1_yp,
            sampler1_x * sampler1_x,
            sampler1_x * sampler1_xp,
            sampler1_x * sampler1_y,
            sampler1_x * sampler1_yp,
            sampler1_xp * sampler1_xp,
            sampler1_xp * sampler1_y,
            sampler1_xp * sampler1_yp,
            sampler1_y * sampler1_y,
            sampler1_y * sampler1_yp,
            sampler1_yp * sampler1_yp,
        ]
    )

    sampler2_matrix = np.array([sampler2_x, sampler2_xp, sampler2_y, sampler2_yp])
    sampler1_matrix_inv = np.linalg.pinv(sampler1_matrix)
    return np.dot(sampler2_matrix, sampler1_matrix_inv)


def map_to_tensor(rmat, tmap_bis):
    rmat6 = np.eye(6)
    rmat6[:4, :4] = np.asarray(rmat, dtype=float)

    tensor = np.zeros((6, 6, 6))
    quad_terms = np.asarray(tmap_bis, dtype=float)[:, 4:]
    for i, row in enumerate(quad_terms):
        for coefficient, (j, k) in zip(row, ORDER2_QUAD_PAIRS):
            if j == k:
                tensor[i, j, k] = coefficient
            else:
                tensor[i, j, k] = coefficient / 2.0
                tensor[i, k, j] = coefficient / 2.0
    return rmat6, tensor


def tensor_to_map(rmat, tensor):
    tmap_bis = np.zeros((4, 14))
    tmap_bis[:, :4] = np.asarray(rmat, dtype=float)[:4, :4]
    for i in range(4):
        for column, (j, k) in enumerate(ORDER2_QUAD_PAIRS, start=4):
            if j == k:
                tmap_bis[i, column] = tensor[i, j, k]
            else:
                tmap_bis[i, column] = tensor[i, j, k] + tensor[i, k, j]
    return tmap_bis


def rotate_taylor_map_order2(rmat, tmap_bis, tilt):
    rmat6, tensor = map_to_tensor(rmat, tmap_bis)
    lab_to_magnet = roll_matrix_bdsim(tilt)
    rmat_lab, tensor_lab = rotate_order2_map_lab_to_magnet(
        rmat6, tensor, lab_to_magnet
    )
    return rmat_lab[:4, :4], tensor_to_map(rmat_lab, tensor_lab)

def run_fint_second_order_consistency_test(
    template_name,
    base_name,
    data,
    fint=0.0,
    fintx=0.0,
    hgap=0.06,
    fintK2=0.0,
    fintxK2=0.0,
    order2_precision=None,
    result_prefix="fint_second_order",
):

    linear_precision=data["linear_precision_fint"]
    if order2_precision is None:
        order2_precision= data["order2_precision_fint"]

    data_zero = data.copy()
    data_zero["TILT1"] = 0.0
    data_zero["E1"] = 0.0
    data_zero["E2"] = 0.0
    data_zero["H1"] = 0.0
    data_zero["H2"] = 0.0
    data_zero["FINT"] = 0.0
    data_zero["FINTX"] = 0.0
    data_zero["FINTK2"] = 0.0
    data_zero["FINTXK2"] = 0.0
    data_zero["HGAP"] = hgap
    
    data_fint = data_zero.copy()
    data_fint["FINT"] = fint
    data_fint["FINTX"] = fintx
    data_fint["FINTK2"] = fintK2
    data_fint["FINTXK2"] = fintxK2

    base_name_zero = base_name + "_zero"
    gmad_name_zero = os.path.join(REGRESSION_DIR, base_name_zero + ".gmad")
    root_name_zero = os.path.join(REGRESSION_DIR, base_name_zero + ".root")
    optics_name_zero = os.path.join(REGRESSION_DIR, base_name_zero + "_optics.root")

    base_name_fint = base_name + "_fint"
    gmad_name_fint = os.path.join(REGRESSION_DIR, base_name_fint + ".gmad")
    root_name_fint= os.path.join(REGRESSION_DIR, base_name_fint + ".root")
    optics_name_fint = os.path.join(REGRESSION_DIR, base_name_fint + "_optics.root")

    r_zero, t_zero = run_order2_case(
        template_name,
        root_name_zero,
        gmad_name_zero,
        optics_name_zero,
        base_name_zero,
        data_zero,
    )

    r_ref0, t_ref0 = build_reference_order2_map(
        data_zero
    )

    body_T_diff = t_zero[:, 4:] - t_ref0[:, 4:]
    #print("body-only T max diff =", np.max(np.abs(body_T_diff)))

    r_fint, t_fint = run_order2_case(
        template_name,
        root_name_fint,
        gmad_name_fint,
        optics_name_fint,
        base_name_fint,
        data_fint,
    )

    r_ref_fint, t_ref_fint = build_reference_order2_map(
        data_fint
    )

    delta_R_bdsim = r_fint - r_zero
    delta_R_ref = r_ref_fint - r_ref0

    delta_T_bdsim = t_fint[:, 4:] - t_zero[:, 4:]
    delta_T_ref = t_ref_fint[:, 4:] - t_ref0[:, 4:]

    r_diff = delta_R_bdsim - delta_R_ref
    t_diff = delta_T_bdsim - delta_T_ref

    print("ΔR BDSIM =", flush=True)
    print(round_matrix(delta_R_bdsim, 6), flush=True)

    print("ΔR reference =", flush=True)
    print(round_matrix(delta_R_ref, 6), flush=True)

    print("ΔR diff =", flush=True)
    print(round_matrix(r_diff, 6), flush=True)

    print("ΔT BDSIM =", flush=True)
    print(round_matrix(delta_T_bdsim, 6), flush=True)

    print("ΔT reference =", flush=True)
    print(round_matrix(delta_T_ref, 6), flush=True)

    print("ΔT diff =", flush=True)
    print(round_matrix(t_diff, 6), flush=True)

    max_r_diff = np.max(np.abs(r_diff))
    max_t_diff = np.max(np.abs(t_diff))

    append_simulation_result(
        data,
        f"{result_prefix}_linear",
        max_r_diff,
        linear_precision,
    )
    append_simulation_result(
        data,
        f"{result_prefix}_order2",
        max_t_diff,
        order2_precision,
    )
    print("FINT test: linear diff: {} must be lower than {}; second order diff: {} must be lower than {}".format(max_r_diff, linear_precision, max_t_diff, order2_precision))

    assert max_r_diff <= linear_precision
    assert max_t_diff <= order2_precision


def build_reference_order2_map(data):
    
    length = data["LENGTH"]
    angle = data["ANGLE"]
    tilt = data["TILT1"]
    e1= data["E1"]
    e2 = data["E2"]
    h1 = data["H1"]
    h2 = data["H2"]
    fint = data["FINT"]
    fintx = data["FINTX"]
    fintK2 = data["FINTK2"]
    fintxK2 = data["FINTXK2"]
    hgap = data["HGAP"]
    body_k1 = data["BODY_K1"]

    rho = length / angle
    
    r_body = np.eye(6)
    r_body[:4, :4] = np.array(
        [
            [np.cos(angle), rho * np.sin(angle), 0.0, 0.0],
            [-(1.0 / rho) * np.sin(angle), np.cos(angle), 0.0, 0.0],
            [0.0, 0.0, 1.0, length],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    
    t_body = sector_body_T_minimal(
        length=length,
        angle=angle,
        body_k1=body_k1,
    )

    r_in = edge_R_matrix(
        rho,
        psi=e1,
        fint=fint,
        hgap=hgap,
    )

    t_in = edge_T_matrix(
        rho,
        psi=e1,
        fint=fint,
        hgap=hgap,
        pole_face_curvature=h1,
        fintK2=fintK2,          
        entrance=True,
    )

    r_out = edge_R_matrix(
        rho,
        psi=e2,
        fint=fintx,
        hgap=hgap,
    )

    t_out = edge_T_matrix(
        rho,
        psi=e2,
        fint=fintx,
        hgap=hgap,
        pole_face_curvature=h2,
        fintK2=fintxK2,        
        entrance=False,
    )

    r_tmp, t_tmp = compose_order2(r_body, t_body, r_in, t_in)
    r_magnet, t_magnet = compose_order2(r_out, t_out, r_tmp, t_tmp)

    lab_to_magnet = roll_matrix_bdsim(tilt)
    r_lab, t_lab = rotate_order2_map_lab_to_magnet(
        r_magnet, t_magnet, lab_to_magnet
    )
    return r_lab[:4, :4], tensor_to_map(r_lab, t_lab)


def run_h1_second_order_difference_test(template_name, base_name, data, h1=0.0, h2=0.0):
    order2_precision=data["order2_precision_h1"]
    data_h1_zero = data.copy()
    data_h1_zero["H1"] = 0.0
    data_h1_zero["H2"] = 0.0
    data_h1_zero["FINT"] = 0.0
    data_h1_zero["FINTX"] = 0.0
    data_h1_zero["FINTK2"] = 0.0
    data_h1_zero["FINTXK2"] = 0.0
    data_h1_zero["TILT1"] = 0.0
    data_h1_nonzero = data.copy()
    data_h1_nonzero["H1"] = h1
    data_h1_nonzero["H2"] = h2
    data_h1_nonzero["TILT1"] = 0.0
    data_h1_nonzero["FINT"] = 0.0
    data_h1_nonzero["FINTX"] = 0.0
    data_h1_nonzero["FINTK2"] = 0.0
    data_h1_nonzero["FINTXK2"] = 0.0

    base_name_zero= base_name + "_h1_zero"
    gmad_name_zero = os.path.join(REGRESSION_DIR, base_name_zero + ".gmad")
    root_name_zero = os.path.join(REGRESSION_DIR, base_name_zero + ".root")
    optics_name_zero = os.path.join(REGRESSION_DIR, base_name_zero + "_optics.root")

    base_name_nonzero= base_name + "_h1_nonzero"
    gmad_name_nonzero = os.path.join(REGRESSION_DIR, base_name_nonzero + ".gmad")
    root_name_nonzero = os.path.join(REGRESSION_DIR, base_name_nonzero + ".root")
    optics_name_nonzero = os.path.join(REGRESSION_DIR, base_name_nonzero + "_optics.root")

    r_zero, t_zero = run_order2_case(template_name,root_name_zero,gmad_name_zero,optics_name_zero,base_name_zero, data_h1_zero)
    r_h1, t_h1 = run_order2_case(template_name,root_name_nonzero,gmad_name_nonzero,optics_name_nonzero,base_name_nonzero, data_h1_nonzero)

    delta_T_bdsim = t_h1[:, 4:] - t_zero[:, 4:]

    _, ref_tmap_zero = build_reference_order2_map(
    data_h1_zero
    )

    _, ref_tmap_h1 = build_reference_order2_map(
        data_h1_nonzero
    )

    delta_T_ref = ref_tmap_h1[:, 4:] - ref_tmap_zero[:, 4:]

    print("ΔR from H1 max bdsim =", np.max(np.abs(r_h1 - r_zero)), flush=True)
    print("ΔT from H1 max bdsim =", np.max(np.abs(delta_T_bdsim)), flush=True)

    print("ΔT from H1 max ref =", np.max(np.abs(delta_T_ref)), flush=True)
    diff = delta_T_bdsim - delta_T_ref
    diff_max = np.max(np.abs(diff))
    print("H1 second order test: {} must be lower than {}".format(diff_max, order2_precision))
    append_simulation_result(data, "h1_second_order", diff_max, order2_precision)
    assert diff_max <= order2_precision


def run_order2_case(
    template_name,
    root_name,
    gmad_name,
    optics_name,
    base_name,
    data
):
    template_dir = os.path.dirname(template_name)
    template_file = os.path.basename(template_name)
    pybdsim.Run.RenderGmadJinjaTemplate(template_file, gmad_name, data, path=template_dir)
    pybdsim.Run.Bdsim(gmad_name, base_name, 3000, 1, silent=True)
    pybdsim.Run.RebdsimOptics(root_name, optics_name, silent=True)

    rmat = pybdsim.Analysis.CalculateRMatrix(root_name, "d1", "t1", size=4, average=True)
    tmap_bis = CalculateTaylorMapOrder2(
        root_name, "d1", "t1", average=True
    )
    return rmat, tmap_bis


def run_second_order_tilt_consistency_test(template_name, base_name, data):

    tilt = data["TILT1"]
    linear_precision = data["linear_precision_tilt"]
    order2_precision = data["order2_precision"]
    data_tilted = data.copy()
    data_tilted["E1"] = 0.0
    data_tilted["E2"] = 0.0
    data_tilted["H1"] = 0.0
    data_tilted["H2"] = 0.0
    data_tilted["FINT"] = 0.0
    data_tilted["FINTX"] = 0.0
    data_tilted["FINTK2"] = 0.0
    data_tilted["FINTXK2"] = 0.0
    data_tilted["HGAP"] = 0.0
    data_tilted["BODY_K1"] = 0.0
    data_untilted = data_tilted.copy()
    data_untilted["TILT1"] = 0

    base_name_untilted = base_name + "_untilted"
    gmad_name_untilted = os.path.join(REGRESSION_DIR, base_name_untilted + ".gmad")
    root_name_untilted = os.path.join(REGRESSION_DIR, base_name_untilted + ".root")
    optics_name_untilted = os.path.join(REGRESSION_DIR, base_name_untilted + "_optics.root")

    base_name_tilted = base_name + "_tilted"
    gmad_name_tilted = os.path.join(REGRESSION_DIR, base_name_tilted + ".gmad")
    root_name_tilted = os.path.join(REGRESSION_DIR, base_name_tilted + ".root")
    optics_name_tilted = os.path.join(REGRESSION_DIR, base_name_tilted + "_optics.root")

    rmat_untilted, tmap_untilted = run_order2_case(
        template_name,
        root_name_untilted,
        gmad_name_untilted,
        optics_name_untilted,
        base_name_untilted,
        data_untilted
    )
    rmat_tilted, tmap_tilted = run_order2_case(
        template_name,
        root_name_tilted,
        gmad_name_tilted,
        optics_name_tilted,
        base_name_tilted,
        data_tilted
    )
    
    rotated_rmat, rotated_tmap = rotate_taylor_map_order2(
        rmat_untilted, tmap_untilted, tilt
    )

    linear_diff = np.max(np.abs(rmat_tilted - rotated_rmat))
    order2_diff = np.max(np.abs(tmap_tilted[:, 4:] - rotated_tmap[:, 4:]))
    append_simulation_result(
        data,
        "second_order_tilt_linear",
        linear_diff,
        linear_precision,
    )
    append_simulation_result(
        data,
        "second_order_tilt_order2",
        order2_diff,
        order2_precision,
    )
    print("Tilt test: linear diff: {} must be lower than {}; second order diff: {} must be lower than {}".format(linear_diff, linear_precision, order2_diff, order2_precision))

    assert linear_diff <= linear_precision
    assert order2_diff <= order2_precision


def run_first_order_test(template_name, base_name, data):
    data1 = data.copy()
    length = data1["LENGTH"]
    angle = data1["ANGLE"]
    tilt1 = data1["TILT1"]
    e1= data1["E1"]
    e2 = data1["E2"]
    data1["BODY_K1"] = 0.0 
    data1["FINT"] = 0.0
    data1["FINTX"] = 0.0
    data1["FINTK2"] = 0.0
    data1["FINTXK2"] = 0.0
    data1["H1"] = 0.0
    data1["H2"] = 0.0
    linear_precision= data1["linear_precision"]    
    rho = length / angle

    base_name = base_name + "_1stOrder"
    gmad_name = os.path.join(REGRESSION_DIR, base_name + ".gmad")
    root_name = os.path.join(REGRESSION_DIR, base_name + ".root")
    optics_name = os.path.join(REGRESSION_DIR, base_name + "_optics.root")

    template_dir = os.path.dirname(template_name)
    template_file = os.path.basename(template_name)
    pybdsim.Run.RenderGmadJinjaTemplate(template_file, gmad_name, data1, path=template_dir)
    pybdsim.Run.Bdsim(gmad_name, base_name, 3000, 1, silent=True)
    pybdsim.Run.RebdsimOptics(root_name, optics_name, silent=True)
    root_file = pybdsim.Data.Load(root_name)

    rmat = pybdsim.Analysis.CalculateRMatrix(root_name, "d1", "t1", size=4, average=True)
    ref_body_rmat = np.array(
        [
            [np.cos(angle), rho * np.sin(angle), 0.0, 0.0],
            [-(1.0 / rho) * np.sin(angle), np.cos(angle), 0.0, 0.0],
            [0.0, 0.0, 1.0, length],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    ref_rmat6 = np.eye(6)
    ref_rmat6[:4, :4] = ref_body_rmat
    ref_rmat = rotated_matrix(rho, ref_rmat6, tilt=tilt1,e1=e1,e2=e2)[:4, :4]
    max_diff = np.max(np.abs(rmat-ref_rmat))
    print("First order test: {} must be lower than {}".format(max_diff, linear_precision))
    append_simulation_result(data, "first_order", max_diff, linear_precision)
    assert max_diff<=linear_precision


def run_first_order_body(template_name, base_name, data):
    data1 = data.copy()

    length = data1["LENGTH"]
    angle = data1["ANGLE"]

    rho = length / angle
    data1["HORIZONTAL_WIDTH"] = min(0.6, rho / 2.0 - 1.0e-3)
    data1["BODY_K1"] = 0.0 
    data1["FINT"] = 0.0
    data1["FINTX"] = 0.0
    data1["FINTK2"] = 0.0
    data1["FINTXK2"] = 0.0
    data1["H1"] = 0.0
    data1["H2"] = 0.0
    data1["TILT1"]=0.0
    data1["E1"]=0.0
    data1["E2"]=0.0

    linear_precision= data1["linear_precision_body"]
    
    base_name = base_name + "_body"
    gmad_name = os.path.join(REGRESSION_DIR, base_name + ".gmad")
    root_name = os.path.join(REGRESSION_DIR, base_name + ".root")
    optics_name = os.path.join(REGRESSION_DIR, base_name + "_optics.root")
    
    template_dir = os.path.dirname(template_name)
    template_file = os.path.basename(template_name)
    
    pybdsim.Run.RenderGmadJinjaTemplate(template_file, gmad_name, data1, path=template_dir)
    pybdsim.Run.Bdsim(gmad_name, base_name, 3000, 1, silent=True)
    pybdsim.Run.RebdsimOptics(root_name, optics_name, silent=True)

    rmat = pybdsim.Analysis.CalculateRMatrix(root_name, "d1", "t1", size=4, average=True)
    ref_rmat = np.array(
        [
            [np.cos(angle), rho * np.sin(angle), 0.0, 0.0],
            [-(1.0 / rho) * np.sin(angle), np.cos(angle), 0.0, 0.0],
            [0.0, 0.0, 1.0, length],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
   
    
    max_diff = np.max(np.abs(rmat-ref_rmat))
    print("First order body test: linear diff: {} must be lower than {}".format(max_diff, linear_precision))
    append_simulation_result(data, "first_order_body", max_diff, linear_precision)
    assert max_diff <= linear_precision

@pytest.mark.parametrize("config", CONFIGS, ids=make_config_id)
def test(config):
    template_name = os.path.join(REGRESSION_DIR, "sbend.tpl")
    base_name = f"sbend_{make_config_id(config)}"

    run_first_order_body(template_name, base_name, config)
    run_first_order_test(template_name, base_name, config)
    run_second_order_tilt_consistency_test(template_name, base_name, config)
    run_h1_second_order_difference_test(template_name, base_name, config, config["H1"], config["H2"])
    run_fint_second_order_consistency_test(
        template_name,
        base_name,
        config,
        fint=0.3,
        fintx=0.3,
        fintK2=0.0,
        fintxK2=0.0,
        result_prefix="fint_second_order",
    )
    run_fint_second_order_consistency_test(
        template_name,
        base_name,
        config,
        fint=0.3,
        fintx=0.3,
        fintK2=0.3,
        fintxK2=0.3,
        order2_precision=config["order2_precision_fint2"],
        result_prefix="fint_k2_second_order",
    )
