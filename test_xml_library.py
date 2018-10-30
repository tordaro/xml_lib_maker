import xml.etree.ElementTree as et
from scipy.constants import pi
from glob import glob

# Define constants
sea_density = 1025.0
steel_density = 7850.0
mass_reduction = (steel_density - sea_density) / steel_density

# Constant values for all materials
pass_alongs = {
    'viscousRollDamping': '0.0',
    'tangentialDragCoeff': '0.0',
    'rayleighStiffness': '0.0',
    'numvelocities': '0',
    'numWaveHeading': '0',
    'massRadius': '0.0',
    'hullnumPoints': '0',
    'dragCoeffz': '1.2',
    'dragCoeffy': '1.2',
    'constructionDamping': '0.0',
    'closeSurfaceNumPoints': '0',
    'LoadType': 'MORRISON'
}

# Collect all xml paths
xml_paths = glob('./*/*.xml')
roots = [et.parse(path).getroot() for path in xml_paths]


def get_matcoeff(root):
    return float(root[0][2].attrib['materialcoefficient'])


def get_description(root):
    return root[0][0].attrib['des']


def get_weightInAir(root):
    return float(root[0][1].attrib['weightInAir'])


def get_area(root):
    return float(root[0][1].attrib['areal'])


def get_dragArea(root):
    return float(root[0][3].attrib['dragArealy'])


def test_get_dragArea():
    '''Vital to check that get_dragArea does its job'''
    truss_tree = et.parse('./test_lib/24 Tuf-3.xml')
    truss_root = truss_tree.getroot()
    exp = 0.024
    obs = get_dragArea(truss_root)
    assert obs == exp


def test_get_area():
    '''Vital to check that get_area() does its job'''
    truss_tree = et.parse('./test_lib/24 Tuf-3.xml')
    truss_root = truss_tree.getroot()
    exp = 0.0004523893421169302
    obs = get_area(truss_root)
    assert obs == exp


def test_get_weightInAir():
    '''Vital to check that get_weightInAir() does its job'''
    truss_tree = et.parse('./test_lib/24 Tuf-3.xml')
    truss_root = truss_tree.getroot()
    exp = 0.38099999999999995
    obs = get_weightInAir(truss_root)
    assert obs == exp


def test_get_description():
    '''Vital to check that get_description() does its job'''
    truss_tree = et.parse('./test_lib/24 Tuf-3.xml')
    truss_root = truss_tree.getroot()
    exp = '24 mm Tufflex 3-slått'
    obs = get_description(truss_root)
    assert obs == exp


def test_truss_matcoeff():
    '''Vital to check that get_matcoeff() does its job'''
    truss_tree = et.parse('./test_lib/24 Tuf-3.xml')
    truss_root = truss_tree.getroot()
    exp = 3.0
    obs = get_matcoeff(truss_root)
    assert exp == obs


def test_chain_matcoeff():
    '''Vital to check that get_matcoeff() does its job'''
    chain_tree = et.parse('./test_lib/32 AnKj.xml')
    chain_root = chain_tree.getroot()
    exp = 2.0
    obs = get_matcoeff(chain_root)
    assert exp == obs


def test_all_matcoeff():
    '''Check that all materialcoefficients are either 2.0 or 3.0.
    Other tests depend on it.'''
    for root in roots:
        exp1 = 2.0
        exp2 = 3.0
        obs = get_matcoeff(root)
        assert obs == exp1 or obs == exp2


def test_emodulus():
    for root in roots:
        des = get_description(root)
        matcoeff = get_matcoeff(root)
        obs = float(root[0][1].attrib['young'])
        if 'Nylon' in des:
            exp = 5.0e8
            assert obs == exp
        elif matcoeff == 3.0:
            exp = 2.0e9
            assert obs == exp
        elif matcoeff == 2.0:
            exp = 1.1e11
            assert obs == exp


def test_weightinwateroverwritten():
    for root in roots:
        exp = 'true'
        obs = root[0][1].attrib['weightinwateroverwritten']
        assert obs == exp


def test_weightInWater():
    for root in roots:
        matcoeff = get_matcoeff(root)
        obs = float(root[0][1].attrib['weightInWater'])
        if matcoeff == 3.0:
            exp = 0.001
            assert obs == exp
        elif matcoeff == 2.0:
            exp = float(root[0][1].attrib['weightInAir']) * mass_reduction
            assert obs == exp


def test_volumeoverwritten():
    for root in roots:
        exp = 'true'
        obs = root[0][1].attrib['volumeoverwritten']
        assert obs == exp


def test_pretension():
    for root in roots:
        exp = 0.0
        obs = float(root[0][1].attrib['pretension'])
        assert obs == exp


def test_volume():
    for root in roots:
        matcoeff = get_matcoeff(root)
        m_l = get_weightInAir(root)
        area = get_area(root)
        obs = float(root[0][1].attrib['volume'])
        if matcoeff == 3.0:
            exp = area
            assert obs == exp
        elif matcoeff == 2.0:
            exp = m_l / steel_density
            assert obs == exp


def test_noCompressionForces():
    for root in roots:
        exp = 'false'
        obs = root[0][1].attrib['noCompressionForces']
        assert obs == exp


def test_dragArea():
    for root in roots:
        exp = get_dragArea(root)
        obs = float(root[0][3].attrib['dragArealyz'])
        assert obs == exp


def test_area():
    for root in roots:
        matcoeff = get_matcoeff(root)
        diameter = get_dragArea(root)
        obs = float(root[0][1].attrib['areal'])
        if matcoeff == 3.0:
            exp = pi * diameter ** 2 / 4
            assert obs == exp
        elif matcoeff == 2.0:
            exp = 2 * pi * diameter ** 2 / 4
            assert obs == exp


def test_massdensity():
    for root in roots:
        m_l = get_weightInAir(root)
        area = get_area(root)
        exp = m_l / area
        obs = float(root[0][1].attrib['massdensity'])
        assert obs == exp


def test_addedmasscoefflocalz():
    for root in roots:
        matcoeff = get_matcoeff(root)
        obs = float(root[0][1].attrib['addedmasscoefflocalz'])
        if matcoeff == 3.0:
            exp = 0.0
            assert obs == exp
        elif matcoeff == 2.0:
            exp = 1.0
            assert obs == exp


def test_addedmasscoefflocaly():
    for root in roots:
        matcoeff = get_matcoeff(root)
        obs = float(root[0][1].attrib['addedmasscoefflocaly'])
        if matcoeff == 3.0:
            exp = 0.0
            assert obs == exp
        elif matcoeff == 2.0:
            exp = 1.0
            assert obs == exp


def test_trusstype():
    for root in roots:
        exp = 3
        obs = int(root[0][2].attrib['trusstype'])
        assert obs == exp


def test_loadmodel():
    for root in roots:
        for key, val in pass_alongs.items():
            exp = val
            obs = root[0][3].attrib[key]
            assert obs == exp


def test_name():
    for root in roots:
        exp = get_dragArea(root) * 1e3
        obs = float(root[0].attrib['name'].split()[0])
        assert obs == exp


def test_description():
    for root in roots:
        exp = get_dragArea(root) * 1e3
        obs = float(root[0][0].attrib['des'].split()[0])
        assert obs == exp
