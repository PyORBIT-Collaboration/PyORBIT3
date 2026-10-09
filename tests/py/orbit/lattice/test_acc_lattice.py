import copy

import pytest

import orbit.py_linac.lattice.LinacAccNodes as linac_nodes
from orbit.lattice import AccLattice, AccNode
from orbit.matrix_lattice import MATRIX_Lattice
from orbit.py_linac.lattice import Bend, Drift, LinacAccLattice, Sequence
from orbit.teapot import MultipoleTEAPOT, QuadTEAPOT, TEAPOT_Lattice


class _CountingNode(AccNode):
    def __init__(self, name="counting"):
        self.initialize_calls = 0
        super().__init__(name)

    def initialize(self):
        self.initialize_calls += 1


class _RequiredArgsLattice(AccLattice):
    __slots__ = ("slot_state",)

    new_calls = 0
    init_calls = 0

    def __new__(cls, required):
        cls.new_calls += 1
        return super().__new__(cls)

    def __init__(self, required):
        type(self).init_calls += 1
        super().__init__("required-args")
        self.required = required
        self.slot_state = ["slot-state"]
        self.__children = ["subclass-state"]


class _NonCopyable:
    def __deepcopy__(self, memo):
        raise TypeError("payload is not copyable")


class _WaveformRead(Exception):
    pass


class _Waveform:
    def __init__(self, owner):
        self.owner = owner

    def getStrength(self):
        raise _WaveformRead(self.owner)


class _Dispatcher:
    def __init__(self):
        self.state = []


class _SlottedDrift(Drift):
    __slots__ = ("slot_state",)

    def __init__(self, name="drift"):
        super().__init__(name)
        self.slot_state = []


class _FakeBunch:
    def charge(self):
        return 1.0


class _FakeTPB:
    def __init__(self):
        self.fringe_strengths = []

    def bendfringeIN(self, bunch, rho):
        pass

    def bendfringeOUT(self, bunch, rho):
        pass

    def multpfringeIN(self, bunch, pole, strength, skew):
        self.fringe_strengths.append(("in", strength))

    def multpfringeOUT(self, bunch, pole, strength, skew):
        self.fringe_strengths.append(("out", strength))


def _closure_values(function):
    return tuple(cell.cell_contents for cell in (function.__closure__ or ()))


def test_shallow_copy_detaches_lattice_containers_but_shares_nodes():
    lattice = AccLattice("source")
    node = AccNode("node")
    node.setLength(1.25)
    lattice.addNode(node)
    lattice.initialize()

    copied = copy.copy(lattice)

    assert type(copied) is type(lattice)
    assert copied is not lattice
    assert copied.getName() == lattice.getName()
    assert copied.getType() == lattice.getType()
    assert copied.getLength() == lattice.getLength()
    assert copied.isInitialized() is lattice.isInitialized()
    assert copied.getNodes() is not lattice.getNodes()
    assert copied.getNodePositionsDict() is not lattice.getNodePositionsDict()
    assert copied.getNodePositionsDict() == lattice.getNodePositionsDict()
    assert copied.getNodes()[0] is node

    node.setParam("shared", [1])
    assert copied.getNodes()[0].getParam("shared") is node.getParam("shared")

    copied.addNode(AccNode("copy-only"))
    assert [item.getName() for item in lattice.getNodes()] == ["node"]
    assert lattice.isInitialized()

    copied.getNodes().remove(node)
    assert lattice.getNodes() == [node]


def test_deepcopy_clones_graph_and_preserves_aliases_and_cycles():
    lattice = AccLattice("source")
    parent = AccNode("parent")
    child = AccNode("child")
    parent.setLength(1.5)
    parent.setParam("payload", {"values": [1]})
    parent.setParam("lattice", lattice)
    child.setParamsDict(parent.getParamsDict())
    parent.addChildNode(child, AccNode.ENTRANCE)
    lattice.addNode(parent)
    lattice.initialize()

    copied = copy.deepcopy(lattice)
    copied_parent = copied.getNodes()[0]
    copied_child = copied_parent.getChildNodes(AccNode.ENTRANCE)[0]

    assert copied_parent is not parent
    assert copied_child is not child
    assert copied_parent.getParamsDict() is copied_child.getParamsDict()
    assert copied_parent.getParamsDict() is not parent.getParamsDict()
    assert copied_parent.getParam("lattice") is copied

    copied_parent.getParam("payload")["values"].append(2)
    assert parent.getParam("payload") == {"values": [1]}

    copied_position_node = next(iter(copied.getNodePositionsDict()))
    assert copied_position_node is copied_parent
    assert copied.getNodePositionsDict()[copied_parent] == (0.0, 1.5)


def test_copy_protocol_bypasses_subclass_constructors_and_base_name_mangling():
    required = ["state"]
    lattice = _RequiredArgsLattice(required)
    lattice.addNode(AccNode("node"))
    calls_before_copy = (_RequiredArgsLattice.new_calls, _RequiredArgsLattice.init_calls)

    shallow = copy.copy(lattice)
    deep = copy.deepcopy(lattice)

    assert (_RequiredArgsLattice.new_calls, _RequiredArgsLattice.init_calls) == calls_before_copy
    assert type(shallow) is _RequiredArgsLattice
    assert type(deep) is _RequiredArgsLattice
    assert shallow.getNodes() is not lattice.getNodes()
    assert shallow.getNodes()[0] is lattice.getNodes()[0]
    assert deep.getNodes()[0] is not lattice.getNodes()[0]
    assert shallow.required is required
    assert deep.required is not required
    assert shallow.slot_state is lattice.slot_state
    assert deep.slot_state is not lattice.slot_state
    assert shallow._RequiredArgsLattice__children is lattice._RequiredArgsLattice__children
    assert deep._RequiredArgsLattice__children is not lattice._RequiredArgsLattice__children


@pytest.mark.parametrize("initialized", [False, True], ids=["uninitialized", "initialized"])
@pytest.mark.parametrize("copy_operation", [copy.copy, copy.deepcopy], ids=["shallow", "deep"])
def test_copy_preserves_initialization_state_without_initializing(copy_operation, initialized):
    lattice = AccLattice("source")
    node = _CountingNode()
    lattice.addNode(node)
    if initialized:
        lattice.initialize()
    initialize_calls = node.initialize_calls

    copied = copy_operation(lattice)

    assert node.initialize_calls == initialize_calls
    assert copied.isInitialized() is initialized
    assert copied.getLength() == lattice.getLength()
    assert list(copied.getNodePositionsDict().values()) == list(lattice.getNodePositionsDict().values())
    assert copied.getNodes()[0].initialize_calls == initialize_calls


def test_deepcopy_reports_noncopyable_top_level_attribute():
    lattice = AccLattice("source")
    lattice.payload = _NonCopyable()

    with pytest.raises(TypeError, match=r"AccLattice.*payload") as exc_info:
        copy.deepcopy(lattice)

    assert isinstance(exc_info.value.__cause__, TypeError)
    assert str(exc_info.value.__cause__) == "payload is not copyable"


def test_matrix_lattice_deepcopy_fails_with_native_attribute_context():
    lattice = MATRIX_Lattice("matrix")

    with pytest.raises(TypeError, match=r"MATRIX_Lattice.*oneTurnMatrix") as exc_info:
        copy.deepcopy(lattice)

    assert isinstance(exc_info.value.__cause__, TypeError)


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(
            lambda waveform: MultipoleTEAPOT(
                "multipole",
                length=1.0,
                poles=[2],
                kls=[1.0],
                skews=[0],
                waveform=waveform,
            ),
            id="multipole",
        ),
        pytest.param(lambda waveform: QuadTEAPOT("quad", length=1.0, waveform=waveform), id="quad"),
    ],
)
@pytest.mark.parametrize(
    ("callback_getter", "fringe_getter"),
    [
        pytest.param("getFringeFieldFunctionIN", "getNodeFringeFieldIN", id="entrance"),
        pytest.param("getFringeFieldFunctionOUT", "getNodeFringeFieldOUT", id="exit"),
    ],
)
def test_deepcopy_teapot_fringe_uses_copied_parent(factory, callback_getter, fringe_getter):
    lattice = TEAPOT_Lattice("teapot")
    source_node = factory(_Waveform("source"))
    lattice.addNode(source_node)
    lattice.initialize()

    copied_lattice = copy.deepcopy(lattice)
    copied_node = copied_lattice.getNodes()[0]
    copied_fringe = getattr(copied_node, fringe_getter)()
    callback = getattr(copied_node, callback_getter)()
    copied_node.waveform.owner = "copy"

    assert copied_node.getParamsDict() is copied_fringe.getParamsDict()
    assert copied_node.getParamsDict() is not source_node.getParamsDict()
    assert source_node not in _closure_values(callback)
    assert next(iter(copied_lattice.getNodePositionsDict())) is copied_node

    with pytest.raises(_WaveformRead, match="copy"):
        callback(copied_fringe, {"parentNode": copied_node, "bunch": _FakeBunch()})


def test_deepcopy_linac_graph_and_bend_callbacks_use_copied_state(monkeypatch):
    lattice = LinacAccLattice("linac")
    sequence = Sequence("sequence")
    sequence.setLinacAccLattice(lattice)
    source_bend = Bend("bend")
    source_bend.setLength(2.0)
    source_bend.setParam("poles", [2])
    source_bend.setParam("kls", [3.0])
    source_bend.setParam("skews", [0])
    sequence.addNode(source_bend)
    lattice.addNode(source_bend)
    lattice.initialize()

    copied_lattice = copy.deepcopy(lattice)
    copied_bend = copied_lattice.getNodes()[0]
    copied_sequence = copied_lattice.getSequences()[0]

    assert copied_bend is not source_bend
    assert copied_sequence is not sequence
    assert copied_bend.getSequence() is copied_sequence
    assert copied_sequence.getNodes()[0] is copied_bend
    assert copied_sequence.getLinacAccLattice() is copied_lattice
    assert copied_bend.tracking_module is source_bend.tracking_module
    assert next(iter(copied_lattice.getNodePositionsDict())) is copied_bend

    copied_bend.setParam("kls", [7.0])
    fake_tpb = _FakeTPB()
    monkeypatch.setattr(linac_nodes, "TPB", fake_tpb)
    params = {"bunch": _FakeBunch(), "parentNode": copied_bend}

    entrance_callback = copied_bend.getFringeFieldFunctionIN()
    exit_callback = copied_bend.getFringeFieldFunctionOUT()
    assert source_bend not in _closure_values(entrance_callback)
    assert source_bend not in _closure_values(exit_callback)
    assert copied_bend.getParamsDict() is copied_bend.getNodeFringeFieldIN().getParamsDict()
    assert copied_bend.getParamsDict() is copied_bend.getNodeFringeFieldOUT().getParamsDict()

    entrance_callback(copied_bend.getNodeFringeFieldIN(), params)
    exit_callback(copied_bend.getNodeFringeFieldOUT(), params)

    assert fake_tpb.fringe_strengths == [("in", -7.0), ("out", -7.0)]


def test_deepcopy_linac_node_recursively_copies_custom_dispatcher():
    node = _SlottedDrift("drift")
    node.tracking_module = _Dispatcher()

    copied = copy.deepcopy(node)

    assert copied.tracking_module is not node.tracking_module
    assert copied.tracking_module.state is not node.tracking_module.state
    assert copied.slot_state is not node.slot_state
