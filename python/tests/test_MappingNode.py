import numpy as np
import pytest

import libcasm.mapping.info as mapinfo
import libcasm.mapping.mapsearch as mapsearch
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims


@pytest.fixture
def hcp_prim():
    """Fixture providing an HCP prim for testing."""
    disp_dof = xtal.DoFSetBasis("disp")
    Hstrain_dof = xtal.DoFSetBasis("Hstrain")
    return xtal_prims.HCP(
        a=1.0,
        occ_dof=["A"],
        local_dof=[disp_dof],
        global_dof=[Hstrain_dof],
    )


@pytest.fixture
def hcp_search_data(hcp_prim):
    """Fixture providing search data objects for testing."""
    hcp_prim_data = mapsearch.PrimSearchData(
        prim=hcp_prim,
        enable_symmetry_breaking_atom_cost=True,
    )
    hcp_prim_structure_data = mapsearch.StructureSearchData(
        lattice=hcp_prim.lattice(),
        atom_coordinate_cart=hcp_prim.coordinate_cart(),
        atom_type=[occ[0] for occ in hcp_prim.occ_dof()],
        override_structure_factor_group=None,
    )
    return hcp_prim_data, hcp_prim_structure_data


@pytest.fixture
def mapping_search():
    """Fixture providing a basic MappingSearch object."""
    return mapsearch.MappingSearch(
        min_cost=0.0,
        max_cost=1e20,
        k_best=10,
        atom_cost_f=mapsearch.IsotropicAtomCost(),
        total_cost_f=mapsearch.WeightedTotalCost(lattice_cost_weight=0.5),
        atom_to_site_cost_f=mapsearch.make_atom_to_site_cost,
        enable_remove_mean_displacement=False,
        infinity=1e20,
        cost_tol=1e-5,
    )


@pytest.fixture
def hcp_mapping_node(hcp_search_data, mapping_search):
    """Fixture providing a basic MappingNode for testing."""
    hcp_prim_data, hcp_prim_structure_data = hcp_search_data
    transformation_matrix_to_super = np.array(
        [
            [1, 0, 0],
            [1, 2, 0],
            [0, 0, 1],
        ],
        dtype="int",
    )

    child_search_data = mapsearch.make_superstructure_data(
        prim_structure_data=hcp_prim_structure_data,
        transformation_matrix_to_super=transformation_matrix_to_super,
    )

    lattice_mappings = mapmethods.map_lattices(
        lattice1=hcp_prim_data.prim_lattice(),
        lattice2=child_search_data.lattice(),
        transformation_matrix_to_super=child_search_data.transformation_matrix_to_super(),
        lattice1_point_group=hcp_prim_data.prim_crystal_point_group(),
        lattice2_point_group=child_search_data.structure_crystal_point_group(),
        min_cost=0.0,
        max_cost=1e20,
        cost_method="isotropic_strain_cost",
        k_best=100,
        reorientation_range=3,
        cost_tol=1e-5,
    )

    scored_lattice_mapping = lattice_mappings[0]
    lattice_mapping_data = mapsearch.LatticeMappingSearchData(
        prim_data=hcp_prim_data,
        structure_data=child_search_data,
        lattice_mapping=scored_lattice_mapping,
    )

    trial_translations = mapsearch.make_trial_translations(
        lattice_mapping_data=lattice_mapping_data,
    )

    # Insert first mapping node
    mapping_search.make_and_insert_mapping_node(
        lattice_cost=scored_lattice_mapping.lattice_cost(),
        lattice_mapping_data=lattice_mapping_data,
        trial_translation_cart=trial_translations[0],
        forced_on={},
        forced_off=[],
    )

    return mapping_search.front()


def test_mapping_node_constructor(hcp_search_data, mapping_search):
    """Test MappingNode constructor."""
    hcp_prim_data, hcp_prim_structure_data = hcp_search_data
    transformation_matrix_to_super = np.array(
        [[1, 0, 0], [1, 2, 0], [0, 0, 1]], dtype="int"
    )

    child_search_data = mapsearch.make_superstructure_data(
        prim_structure_data=hcp_prim_structure_data,
        transformation_matrix_to_super=transformation_matrix_to_super,
    )

    lattice_mappings = mapmethods.map_lattices(
        lattice1=hcp_prim_data.prim_lattice(),
        lattice2=child_search_data.lattice(),
        transformation_matrix_to_super=child_search_data.transformation_matrix_to_super(),
        lattice1_point_group=hcp_prim_data.prim_crystal_point_group(),
        lattice2_point_group=child_search_data.structure_crystal_point_group(),
        min_cost=0.0,
        max_cost=1e20,
        cost_method="isotropic_strain_cost",
        k_best=1,
        reorientation_range=1,
        cost_tol=1e-5,
    )

    lattice_mapping_data = mapsearch.LatticeMappingSearchData(
        prim_data=hcp_prim_data,
        structure_data=child_search_data,
        lattice_mapping=lattice_mappings[0],
    )

    trial_translation = np.array([0.0, 0.0, 0.0])

    # Test constructor with default forced_on and forced_off
    node = mapsearch.MappingNode(
        search=mapping_search,
        lattice_cost=lattice_mappings[0].lattice_cost(),
        lattice_mapping_data=lattice_mapping_data,
        trial_translation_cart=trial_translation,
        forced_on={},
        forced_off=[],
    )

    assert node is not None
    assert isinstance(node.lattice_cost(), float)
    assert isinstance(node.atom_cost(), float)
    assert isinstance(node.total_cost(), float)


def test_mapping_node_lattice_cost(hcp_mapping_node):
    """Test MappingNode.lattice_cost method."""
    assert hasattr(hcp_mapping_node, "lattice_cost")
    lattice_cost = hcp_mapping_node.lattice_cost()
    assert isinstance(lattice_cost, float)
    assert lattice_cost >= 0.0


def test_mapping_node_atom_cost(hcp_mapping_node):
    """Test MappingNode.atom_cost method."""
    assert hasattr(hcp_mapping_node, "atom_cost")
    atom_cost = hcp_mapping_node.atom_cost()
    assert isinstance(atom_cost, float)
    assert atom_cost >= 0.0


def test_mapping_node_total_cost(hcp_mapping_node):
    """Test MappingNode.total_cost method."""
    assert hasattr(hcp_mapping_node, "total_cost")
    total_cost = hcp_mapping_node.total_cost()
    assert isinstance(total_cost, float)
    assert total_cost >= 0.0


def test_mapping_node_lattice_mapping_data(hcp_mapping_node):
    """Test MappingNode.lattice_mapping_data attribute."""
    assert hasattr(hcp_mapping_node, "lattice_mapping_data")
    lattice_mapping_data = hcp_mapping_node.lattice_mapping_data()
    assert isinstance(lattice_mapping_data, mapsearch.LatticeMappingSearchData)


def test_mapping_node_atom_mapping_data(hcp_mapping_node):
    """Test MappingNode.atom_mapping_data attribute."""
    assert hasattr(hcp_mapping_node, "atom_mapping_data")
    atom_mapping_data = hcp_mapping_node.atom_mapping_data()
    assert isinstance(atom_mapping_data, mapsearch.AtomMappingSearchData)


def test_mapping_node_atom_mapping(hcp_mapping_node):
    """Test MappingNode.atom_mapping attribute."""
    assert hasattr(hcp_mapping_node, "atom_mapping")
    atom_mapping = hcp_mapping_node.atom_mapping()
    assert isinstance(atom_mapping, mapinfo.AtomMapping)


def test_mapping_node_forced_on(hcp_mapping_node):
    """Test MappingNode.forced_on attribute."""
    assert hasattr(hcp_mapping_node, "forced_on")
    forced_on = hcp_mapping_node.forced_on()
    assert isinstance(forced_on, dict)


def test_mapping_node_forced_off(hcp_mapping_node):
    """Test MappingNode.forced_off attribute."""
    assert hasattr(hcp_mapping_node, "forced_off")
    forced_off = hcp_mapping_node.forced_off()
    assert isinstance(forced_off, list)


def test_mapping_node_assignment(hcp_mapping_node):
    """Test MappingNode.assignment method."""
    assert hasattr(hcp_mapping_node, "assignment")
    assignment = hcp_mapping_node.assignment()
    assert isinstance(assignment, list)
    assert len(assignment) > 0
    # Check that all assignments are valid indices
    for idx in assignment:
        assert isinstance(idx, int)
        assert idx >= 0


def test_mapping_node_assignment_cost(hcp_mapping_node):
    """Test MappingNode.assignment_cost attribute."""
    assert hasattr(hcp_mapping_node, "assignment_cost")
    assignment_cost = hcp_mapping_node.assignment_cost()
    assert isinstance(assignment_cost, float)
    assert assignment_cost >= 0.0


def test_mapping_node_symmetry_preserving_displacement(hcp_mapping_node):
    """Test MappingNode.symmetry_preserving_displacement method."""
    assert hasattr(hcp_mapping_node, "symmetry_preserving_displacement")
    disp = hcp_mapping_node.symmetry_preserving_displacement()
    assert isinstance(disp, np.ndarray)
    assert disp.shape[0] == 3  # 3D displacement


def test_mapping_node_symmetry_breaking_displacement(hcp_mapping_node):
    """Test MappingNode.symmetry_breaking_displacement method."""
    assert hasattr(hcp_mapping_node, "symmetry_breaking_displacement")
    disp = hcp_mapping_node.symmetry_breaking_displacement()
    assert isinstance(disp, np.ndarray)
    assert disp.shape[0] == 3  # 3D displacement


def test_mapping_node_to_dict(hcp_mapping_node):
    """Test MappingNode.to_dict method."""
    assert hasattr(hcp_mapping_node, "to_dict")
    d = hcp_mapping_node.to_dict()
    assert isinstance(d, dict)

    # Check that required keys are present
    expected_keys = [
        "forced_on",
        "forced_off",
        "atom_mapping",
        "atom_cost",
        "lattice_mapping",
        "lattice_cost",
        "total_cost",
    ]
    for key in expected_keys:
        assert key in d

    # Check types of values
    assert isinstance(d["forced_on"], list)
    assert isinstance(d["forced_off"], list)
    assert isinstance(d["atom_mapping"], dict)
    assert isinstance(d["atom_cost"], float)
    assert isinstance(d["lattice_mapping"], dict)
    assert isinstance(d["lattice_cost"], float)
    assert isinstance(d["total_cost"], float)


def test_mapping_node_with_forced_on(hcp_search_data, mapping_search):
    """Test MappingNode with forced_on constraints."""
    hcp_prim_data, hcp_prim_structure_data = hcp_search_data
    transformation_matrix_to_super = np.array(
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype="int"
    )

    child_search_data = mapsearch.make_superstructure_data(
        prim_structure_data=hcp_prim_structure_data,
        transformation_matrix_to_super=transformation_matrix_to_super,
    )

    lattice_mappings = mapmethods.map_lattices(
        lattice1=hcp_prim_data.prim_lattice(),
        lattice2=child_search_data.lattice(),
        transformation_matrix_to_super=child_search_data.transformation_matrix_to_super(),
        lattice1_point_group=hcp_prim_data.prim_crystal_point_group(),
        lattice2_point_group=child_search_data.structure_crystal_point_group(),
        min_cost=0.0,
        max_cost=1e20,
        cost_method="isotropic_strain_cost",
        k_best=1,
        reorientation_range=1,
        cost_tol=1e-5,
    )

    lattice_mapping_data = mapsearch.LatticeMappingSearchData(
        prim_data=hcp_prim_data,
        structure_data=child_search_data,
        lattice_mapping=lattice_mappings[0],
    )

    trial_translation = np.array([0.0, 0.0, 0.0])
    forced_on = {0: 0}  # Force site 0 to map to atom 0

    node = mapsearch.MappingNode(
        search=mapping_search,
        lattice_cost=lattice_mappings[0].lattice_cost(),
        lattice_mapping_data=lattice_mapping_data,
        trial_translation_cart=trial_translation,
        forced_on=forced_on,
        forced_off=[],
    )

    assert node.forced_on() == forced_on
    assignment = node.assignment()
    assert assignment[0] == 0  # Verify forced assignment


def test_mapping_node_with_forced_off(hcp_search_data, mapping_search):
    """Test MappingNode with forced_off constraints."""
    hcp_prim_data, hcp_prim_structure_data = hcp_search_data
    transformation_matrix_to_super = np.array(
        [[2, 0, 0], [0, 2, 0], [0, 0, 2]], dtype="int"
    )

    child_search_data = mapsearch.make_superstructure_data(
        prim_structure_data=hcp_prim_structure_data,
        transformation_matrix_to_super=transformation_matrix_to_super,
    )

    lattice_mappings = mapmethods.map_lattices(
        lattice1=hcp_prim_data.prim_lattice(),
        lattice2=child_search_data.lattice(),
        transformation_matrix_to_super=child_search_data.transformation_matrix_to_super(),
        lattice1_point_group=hcp_prim_data.prim_crystal_point_group(),
        lattice2_point_group=child_search_data.structure_crystal_point_group(),
        min_cost=0.0,
        max_cost=1e20,
        cost_method="isotropic_strain_cost",
        k_best=1,
        reorientation_range=1,
        cost_tol=1e-5,
    )

    lattice_mapping_data = mapsearch.LatticeMappingSearchData(
        prim_data=hcp_prim_data,
        structure_data=child_search_data,
        lattice_mapping=lattice_mappings[0],
    )

    trial_translation = np.array([0.0, 0.0, 0.0])
    forced_off = [(0, 1)]  # Prevent site 0 from mapping to atom 1

    node = mapsearch.MappingNode(
        search=mapping_search,
        lattice_cost=lattice_mappings[0].lattice_cost(),
        lattice_mapping_data=lattice_mapping_data,
        trial_translation_cart=trial_translation,
        forced_on={},
        forced_off=forced_off,
    )

    assert node.forced_off() == forced_off
    assignment = node.assignment()
    # Verify forced_off constraint: site 0 should not map to atom 1
    if len(assignment) > 0:
        assert assignment[0] != 1


def test_mapping_node_displacement_decomposition(hcp_mapping_node):
    """Test that symmetry-preserving and symmetry-breaking displacements sum
    correctly."""
    atom_mapping = hcp_mapping_node.atom_mapping()
    total_displacement = atom_mapping.displacement()

    sym_preserving = hcp_mapping_node.symmetry_preserving_displacement()
    sym_breaking = hcp_mapping_node.symmetry_breaking_displacement()

    # Check shapes match
    assert sym_preserving.shape == total_displacement.shape
    assert sym_breaking.shape == total_displacement.shape

    # The sum of symmetry-preserving and symmetry-breaking should equal total
    reconstructed = sym_preserving + sym_breaking
    assert np.allclose(reconstructed, total_displacement, atol=1e-10)


def test_mapping_node_comparison(hcp_search_data, mapping_search):
    """Test MappingNode comparison based on total_cost."""
    hcp_prim_data, hcp_prim_structure_data = hcp_search_data
    transformation_matrix_to_super = np.array(
        [[1, 0, 0], [1, 2, 0], [0, 0, 1]], dtype="int"
    )

    child_search_data = mapsearch.make_superstructure_data(
        prim_structure_data=hcp_prim_structure_data,
        transformation_matrix_to_super=transformation_matrix_to_super,
    )

    lattice_mappings = mapmethods.map_lattices(
        lattice1=hcp_prim_data.prim_lattice(),
        lattice2=child_search_data.lattice(),
        transformation_matrix_to_super=child_search_data.transformation_matrix_to_super(),
        lattice1_point_group=hcp_prim_data.prim_crystal_point_group(),
        lattice2_point_group=child_search_data.structure_crystal_point_group(),
        min_cost=0.0,
        max_cost=1e20,
        cost_method="isotropic_strain_cost",
        k_best=10,
        reorientation_range=2,
        cost_tol=1e-5,
    )

    # Create mapping search and insert multiple nodes
    search = mapsearch.MappingSearch(
        min_cost=0.0,
        max_cost=1e20,
        k_best=100,
        atom_cost_f=mapsearch.IsotropicAtomCost(),
        total_cost_f=mapsearch.WeightedTotalCost(lattice_cost_weight=0.5),
        atom_to_site_cost_f=mapsearch.make_atom_to_site_cost,
        enable_remove_mean_displacement=False,
        infinity=1e20,
        cost_tol=1e-5,
    )

    # Iterate through first 3 lattice mappings
    count = 0
    for scored_lattice_mapping in lattice_mappings:
        if count >= 3:
            break
        lattice_mapping_data = mapsearch.LatticeMappingSearchData(
            prim_data=hcp_prim_data,
            structure_data=child_search_data,
            lattice_mapping=scored_lattice_mapping,
        )

        trial_translations = mapsearch.make_trial_translations(
            lattice_mapping_data=lattice_mapping_data,
        )

        search.make_and_insert_mapping_node(
            lattice_cost=scored_lattice_mapping.lattice_cost(),
            lattice_mapping_data=lattice_mapping_data,
            trial_translation_cart=trial_translations[0],
            forced_on={},
            forced_off=[],
        )
        count += 1

    assert search.size() >= 2

    # Get front and back nodes
    front_node = search.front()
    back_node = search.back()

    # Front should have lower or equal cost compared to back
    assert front_node.total_cost() <= back_node.total_cost()


def test_mapping_node_consistency_with_search_data(hcp_mapping_node):
    """Test that MappingNode data is consistent with its search data."""
    lattice_mapping_data = hcp_mapping_node.lattice_mapping_data()
    atom_mapping_data = hcp_mapping_node.atom_mapping_data()
    atom_mapping = hcp_mapping_node.atom_mapping()

    # Check that lattice_mapping_data is consistent
    assert lattice_mapping_data.prim_data() is not None
    assert lattice_mapping_data.structure_data() is not None

    # Check that atom_mapping_data references the same lattice_mapping_data
    assert atom_mapping_data.lattice_mapping_data() is lattice_mapping_data

    # Check that assignment length matches the number of sites
    assignment = hcp_mapping_node.assignment()
    N_sites = lattice_mapping_data.N_supercell_site()
    assert len(assignment) == N_sites

    # Check atom_mapping type consistency
    assert isinstance(atom_mapping, mapinfo.AtomMapping)
