# These classes are meant as containers for general information about pieces of
# secondary structure, as well as utility methods for finding and assessing the
# continued chemical properties of these elements.
import pyrosetta.rosetta as _pyr
from itertools import groupby as _gb
from ..util import or_compose_residue_selectors as or_combine


class SecondaryStructureContainerFactory:
    def __init__(self):
        self._creators = {}

    def register_format(self, dssp_type, creator):
        self._creators[dssp_type] = creator

    def get_container(self, pose, start_pos, end_pos, dssp_type):
        creator = self._creators.get(dssp_type)
        if not creator:
            return SecondaryStructureResidueContainer(
                pose, start_pos, end_pos, dssp_type
            )
        return creator(pose, start_pos, end_pos, dssp_type)


class SecondaryStructureResidueContainer(object):
    """
    An informational container for rosetta pose elements of secondary structure

    Designed to be subclassed into "helix", "strand","gamma turn", etc etc

    This container is not designed to be dynamic and continue to store imformation about the structure over time! Subclasses should have methods to determine continued accuracy of the container for the pose to which it is bound.
    #TODO apply_label()
    #TODO NotImplementedError() for check_residues
    #TODO create_subpose()
    """

    def __init__(self, pose, start_pos, end_pos, dssp_type):
        self.pose = pose
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.residues = pose.residues[start_pos:end_pos]
        self.dssp_type = dssp_type

    def __repr__(self):
        return str(
            {
                "type": type(self),
                "pose": self.pose,
                "start_pos": self.start_pos,
                "end_pos": self.end_pos,
                "dssp_type": self.dssp_type,
            }
        )

    def __str__(self):
        s = (
            f"""{type(self).__name__}:\n"""
            f"""Pose:\n{self.pose}\n"""
            f"""start_pos: {self.start_pos}\n"""
            f"""end_pos: {self.end_pos}\n"""
            f"""dssp_type: {self.dssp_type}"""
        )
        return s

    def get_range(self):
        """
        Returns a tuple (start_pos,end_pos)
        """
        return self.start_pos, self.end_pos

    def resnum_list(self, upstream=0, downstream=0):
        """
        returns a list of residues that are in the container
        """
        return [
            *range(self.start_pos - upstream, self.end_pos + downstream + 1)
        ]

    # TODO maybe get rid of this shift
    def frame_shift(self, n):
        """
        Add n to both the start and end of the SecondaryStructure (accepts neg)
        """
        self.start_pos, self.end_pos = self.start_pos + n, self.end_pos + n

    # TODO: split_structure method, where the structure can be broken at some residue into two identical structures
    def split_structure(self, index):
        """
        return as tuple two strucures: (start,index) & (index +1,end)

        Should error if you try to split somewhere ~unusual~
        """
        return (
            SecondaryStructureResidueContainer(self.start_pos, index),
            SecondaryStructureResidueContainer(index + 1, self.end_pos),
        )

    def generate_residue_selector(self, upstream=0, downstream=0):
        """
        Returns a residue selector for the given container

        upstream and downstream specify a range to allow before and after the
        container. Negative values reduce the selection of the selector.
        """
        start = self.start_pos - upstream
        end = self.end_pos + downstream
        assert bool(
            start <= end
        ), f"Cannot generate residue selector where end: {end} precedes start: {start}"
        if start < 1 or end > len(self.pose.residues):
            raise (
                ValueError,
                f"start: {start} or end: {end} are not found in the pose",
            )

        return _pyr.core.select.residue_selector.ResidueIndexSelector(
            f"{self.start_pos - upstream}-{self.end_pos + downstream}"
        )

    def subpose_structure(self):
        """
        Returns an instance of this where pose only has residues in this container
        """
        subpose = _pyr.protocols.grafting.return_region(
            self.pose.clone(), self.start_pos, self.end_pos
        )
        return type(self)(subpose, 1, len(subpose.residues), self.dssp_type)


class HelixResidueContainer(SecondaryStructureResidueContainer):
    """
    A container to describe consecutive helical residues

    TODO: support for kinks, pitch, coiled coil pairs/indexing, check helicity
    """

    def __init__(self, pose, start_pos, end_pos, dssp_type):
        super().__init__(pose, start_pos, end_pos, dssp_type)

    # TODO check_helix
    """
    Checks if all residues start-end are still helical
    """


class LoopResidueContainer(SecondaryStructureResidueContainer):
    """
    A secondary structure container for holding all types of loops

    TODO preceding and following attributes that are bound to secondary structure up and down stream of the loop
    TODO internal container of all structured turns in the loop
    """

    def __init__(self, pose, start_pos, end_pos, dssp_type):
        super().__init__(pose, start_pos, end_pos, dssp_type)


# TODO class Turn(LoopResidueContainer)


class StrandResidueContainer(SecondaryStructureResidueContainer):
    """
    A container to describe consecutive beta strand residues

    TODO: support for bulge, pitch, neighboring strands, check beta
    """

    def __init__(self, pose, start_pos, end_pos, dssp_type):
        super().__init__(pose, start_pos, end_pos, dssp_type)

    # TODO check_beta
    """
    Checks if all residues start-end are still beta character
    """


def get_container_creator():
    factory = SecondaryStructureContainerFactory()
    factory.register_format("E", StrandResidueContainer)
    factory.register_format("H", HelixResidueContainer)
    factory.register_format("L", LoopResidueContainer)
    return factory


def get_allowed_dssp_values():
    """
    Returns the DSSP string values registered with the current creator
    """
    from copy import deepcopy

    return deepcopy(get_container_creator()._creators)


def parse_structure_from_dssp(pose, *dssp_types):
    """
    returns a list of secondary structure containers for the given pose

    leaving dssp_types blank allows all types. Unsupported dssp values use the
    generic SecondaryStructureResidueContainer
    """

    dssp_dicts = [
        {
            "pose": chain,
            "dssp_string": _pyr.core.scoring.dssp.Dssp(
                chain
            ).get_dssp_secstruct(),
        }
        for chain in pose.split_by_chain()
    ]
    creator = get_container_creator()
    return [
        creator.get_container(
            dssp_dict["pose"], run[0], run[-1], str(res_type_string)
        )
        for dssp_dict in dssp_dicts
        for res_type_string, iterator in _gb(
            enumerate(dssp_dict["dssp_string"], 1), lambda x: x[1]
        )
        for run in ([n for n, p in iterator],)
        if (not dssp_types) or res_type_string in dssp_types
    ]


def remove_ss_from_pose(pose, *dssp_types, in_place=True):
    """
    Must provide dssp_types
    """
    if not dssp_types:
        raise ValueError("At least one dssp type must be provided")
    struct_list = parse_structure_from_dssp(pose, *dssp_types)
    selectors = [ss.generate_residue_selector() for ss in struct_list]
    selector = or_combine(*selectors)
    delete_mover = _pyr.protocols.grafting.simple_movers.DeleteRegionMover()
    delete_mover.set_residue_selector(selector)
    if in_place:
        delete_mover.apply(pose)
        return pose
    else:
        out = pose.clone()
        delete_mover.apply(out)
        return out


def split_ss_to_subpose_containers(pose, *dssp_types):
    """
    Returns a list of containers referencing subposes with the given residues

    Blank dssp_types assumes all
    """
    struct_list = parse_structure_from_dssp(pose, *dssp_types)
    return [ss.subpose_structure() for ss in struct_list]


def split_ss_to_subposes(pose, *dssp_types):
    """
    Returns a list of subposes split by secondary structure

    Blank dssp_types assumes all
    """
    struct_list = parse_structure_from_dssp(pose, *dssp_types)
    return [ss.subpose_structure().pose for ss in struct_list]
