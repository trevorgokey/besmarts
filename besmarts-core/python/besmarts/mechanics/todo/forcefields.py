"""
besmarts.mm.forcefields
"""

from typing import List, Dict

from besmarts.core import trees, hierarchies
from besmarts.core import mm
from besmarts.core import assignments

# notes for how to do the other things
# kept because it has notes for transcoding
"""
    forces between points

    given a SMARTS structure, give the energy and force

    the force is defined over a explicit set of atoms, such as an internal
    coordinate

    Assignment will occur by going over the list of hierarchies and assign each
    parameter based on

    point n cartesian_conformation ; cmd symbol name
    system_parameter x max_assignment int 10 ; cmd symbol name cast value
    hierarchy a atom x y z ; cmd symbol topology symbols
    parameter 1 1 0.0 0.0 0.0 [@BUTANE.1] ; cmd idx lvl x y z SMARTS
    parameter 2 1 0.0 0.0 0.0 [@BUTANE.2]
    parameter 3 1 0.0 0.0 0.0 [@BUTANE.3]
    parameter 4 1 0.0 0.0 0.0 [@BUTANE.4]
    parameter 5 1 0.0 0.0 0.0 [@BUTANE.5]
    parameter 6 1 0.0 0.0 0.0 [@BUTANE.6]
    parameter 7 1 0.0 0.0 0.0 [@BUTANE.7]
    parameter 8 1 0.0 0.0 0.0 [@BUTANE.8]

    force n coulomb_cutoff_switched
    topology_term e partial_charge e

    hierarchy a atom q ; method am1bcc
    parameter 1 1 0.1 [@BUTANE.1]
    parameter 2 1 0.1 [@BUTANE.2]
    parameter 3 1 0.1 [@BUTANE.3]
    parameter 4 1 0.1 [@BUTANE.4]
    parameter 5 1 0.1 [@BUTANE.5]
    parameter 6 1 0.1 [@BUTANE.6]
    parameter 7 1 0.1 [@BUTANE.7]
    parameter 8 1 0.1 [@BUTANE.8]

    @ is a arbitrary name/map that won't get encoded as SMARTS; orthogonal to SMARTS primitives
    $ is a embedded graph that is assumed to be unioned? intersections can be removed
    % is a condensed SMARTS pattern, assumed OR
    ^ is a condensed SMARTS pattern, assumed AND

    An example:
    force n lennard_jones_berkowitz_cutoff_switched
    system_term c cutoff A 10.0
    system_term l switch A 2
    topology_term e epsilon kJ/mol
    topology_term s sigma   A
    topology_term x scale float

    hierarchy a atom epsilon sigma
    parameter 1 1 .01 1.5 [*:1]   ; comment
    parameter 2 1 .02 1.2 [#6:1]  ; comment

    hierarchy b bond scale
    parameter 1 1 1.0 [*:1].[*:2]  ; comment
    parameter 1 1 0.0 [*:1]~[*:2]  ; comment

    hierarchy c pair scale
    parameter 1 1 0.0 [*:1]~[*]~[*:2]  ; comment

    hierarchy d pair scale
    parameter 1 1 0.5 [*:1]~[*]~[*]~[*:2]  ; comment

    extra_points vs fixed_bond_relative_charge
    topology_term d length A

    we might want to apply the transcoding first and pair it with an ff

    we can assign system, mass, bond, angle, torsion, outofplane, nbpairs

    chemical_model system
        procedure system space_box
            system_term x xlength A float 10.0
            system_term y ylength A float 10.0
            system_term z zlength A float 10.0

            system_term a xorigin A float  0.0
            system_term b yorigin A float  0.0
            system_term c zorigin A float  0.0

        model system temperature
            system_term t initial_temperature K float 100.0

        model system pressure
            system_term t initial_pressure K float 100.0
            
        model system volume_box
            system_term t constant_volume nm^3 float 10000

        model system pme_box
            system_term j xgrid A float .1
            system_term k ygrid A float .1
            system_term l zgrid A float .1

    chemical_system
        chemical_model
            procedure
    
    physical_system


    chemical_system bff

        perception rdkit

        bond_to_angle constrained_points
            assign smarts_bond s w11 w12 w13 w21 w22 w23 r
                1 [#9:1]~[*:2] [*:1]~[*:2]~[#0:3] 1 0 0 0 1 0 1.5

        angle_to_atom pack
            assign smarts_bond s w11 w21 w31
                1 [*:1]~[*:2]~[*:3] [@1:1] 1 1 1

        atom_to_graph unpack
            assign smarts s 
                1 [@1:1] [#1:1]-[#8]-[#1:2]

        bond harmonic
            assign smarts_bond k l
                1 [*:1]~[*:2] 0.0 1.5

        angle harmonic
            assign smarts_angle k l
                1 [*:1]~[*:2]~[*:3] 0.0 1.5

        torsion periodic
            assign smarts_torsion k n p
                1 [*:1]~[*:2]~[*:3]~[*:4] 1.0 1 0.0

        outofplane periodic
            assign smarts_outofplane k n p
                1 [*:1]~[*:2](~[*:3])~[*:4] 1.0 1 0.0

        atom coulomb

            parameter combine geometric

            system_term cutoff float A 10.0

            assign smarts_atom q
                parameter 1 [*:1] 0.0

            modify smarts_pair q
                parameter 1 [*:1]~[*:2] 0.0
                parameter 2 [*:1]~[*]~[*:2] 0.0
                parameter 3 [*:1]~[*]~[*]~[*:2] .83
                parameter 4 [*:1]~[*]~[*]~[*]~[*:2] .5

        atom lennard_jones
            parameter combine berkowitz
            assign smarts_atom e s
                parameter 1 [*:1] .5 1.5
            modify smarts_pair e s
                parameter 1 [*:1]~[*:2] 0.0 0.0
                parameter 2 [*:1]~[*]~[*:2] 0.0 0.0
                parameter 3 [*:1]~[*]~[*]~[*:2] .83 1.0
                parameter 4 [*:1]~[*]~[*]~[*]~[*:2] .5 1.0


    chemical_model bff2
        

    chemical_model nbpairs
        model atom coulomb

            system_term d dielectric none float 1.0
            system_term c cutoff A float 10.0

            topology_term q charge e float

            assign smarts_set q
                parameter 1 [*:1] 0.0

            assign extern q
                parameter program sqm
                parameter method am1bcc
                parameter modify no
                parameter conformations rdkit 100
                parameter filter elf10
            
            modify smarts_bond_charge_transfer q
                parameter 1 [*:1]~[*:2] -.1

            modify smarts_angle_charge_transfer q
                parameter 1 [*:1]~[*:2]~[*:3] -.1

        model vdw_berkowitz atom
            topology_term e depth kJ/mol float
            topology_term s sigma A float

            assign smarts e s
                parameter 1 [*:1] 0.0 0.0

            
    chemical_model mass
        model atom mass
            topology_term m mass dalton float
            assign file m
                parameter load besmarts_mass.bp ; from a file
                parameter 1 [#1:1] 12 # comment ; whats in the file

    chemical_model units
        model unit [@1:1]
            assign unit_bond k l
                parameter bond 1 2 1 1.5
                parameter bond 2 3 1 1.5

            assign unit_angle k l
                parameter bond 1 2 3 1 1.5
                parameter bond 2 3 4 1 1.5

    

    chemical_model transode
            
        procedure transcoder angle_to_atom

            parameter mass_weighted true

            topology_term s code none str
            topology_term a weight1 none float
            topology_term b weight2 none float
            topology_term c weight3 none float

            assign smarts s a b c
                parameter u1 0 [!#0:1]~[!#0:2]~[!#0:3] [#0:1] 0 1 0

        model bond_to_angle

            topology_term s code none str
            topology_term a weight1to1 none float
            topology_term b weight1to2 none float
            topology_term c weight1to3 none float
            topology_term e weight2to1 none float
            topology_term f weight2to2 none float
            topology_term g weight2to3 none float
            topology_term r distance   A    float

            assign smarts s a b c
                parameter u1 0 [!@!#0:1]~[!@!#0:2] [*:1][*:2][#0:3] 1 0 1 0 1 0 2.0

        model angle_to_bond

            topology_term s code none str
            topology_term a weight1to1 none float
            topology_term b weight1to2 none float
            topology_term c weight2to1 none float
            topology_term d weight2to2 none float
            topology_term e weight3to1 none float
            topology_term f weight3to2 none float
            topology_term x distance   A    float
            topology_term y distance   A    float

            assign smarts s a b c
                parameter u1 0 [!@!#0:1]~[!@!#0:2] [*:1][*:2][#0:3] 1 0 1 0 1 0 2.0


    chemical_model bond 


        model bond force_harmonic
            topology_term k stiffness kJ/mol/A/A float 1.5
            topology_term l length A float 1.0
            topology_term w wbo A float 
            topology_term a wbo_min A float 1.0
            topology_term b wbo_max A float 2.0

            assign smarts k l # this is a smarts_assignment_hierarchy and k l are the keys
                parameter b1 0 [*:1]~[*:2] 1 1.5  # comment

            assign wbo w
                parameter method sqm
                parameter modify no
                parameter conformations rdkit 100
                parameter filter elf10
                parameter ensemble average
            
            modify wbo_interp k
                need w a b

    chemical_model angle_harmonic

        force angle_harmonic
            topology_term l angle degrees float
            topology_term k stiffness kJ/mol/rad/rad float

            smarts a angle k l
                parameter 1 [*:1]~[*:2]~[*:3] 1 1.5 # comment

    chemical_model vsite_positions

        force bond_constrain
            topology_term l length A float 1.0

            topology_procedure bond l smarts
                parameter 1 [#0:1]~[*:2] 1.0  # comment

        force angle_constrain
            topology_term l angle degrees float
            topology_procedure bond l smarts
                parameter 1 [#0:1]~[*:2]~[*:3] 180.0  # comment

        force torsion_constrain
            topology_term a angle deg float 1.0

            topology_procedure torsion a smarts
                parameter 1 [#0:1]~[*:2]~[*:3]~[*:4] 180  # comment

        force outofplane_constrain
            topology_term a angle deg float 1.0

            topology_procedure outofplane a smarts
                parameter 1 [#0:1]~[*:2](~[*:3])~[*:4] 0  # comment

    chemical_model torsion_periodic

        force torsion_periodic
            topology_term n periodicity none int
            topology_term k amplitude kJ/mol float
            topology_term p phase rad float

            smarts a torsion n k p
                parameter 1 [*:1]~[*:2]~[*:3]~[*:4] 1 1.0 1.0 2 1.0 0.0 # comment

    chemical_model outofplane_periodic

        force outofplane_periodic
            topology_term n periodicity none int
            topology_term k amplitude kJ/mol float
            topology_term p phase rad float

            smarts a outofplane n k p
                parameter 1 [*:1]~[*:2](~[*:3])~[*:4] 1 1 0.1 2 1.0 0.0 # comment

    chemical_model transcodes
        transcoder t angle_to_torsion
            hierarchy a
                parameter -1.3 120 180 [#1:1]-[#8:2]-[#1:3]-[#0:4] ; #0 is a lone pair

    chemical_model elecrostatics

        force coulomb_cutoff
            system_term c cutoff A float 10.0
            system_term e dielectric A float 1.0

            topology_term q charge e float

            topology_procedure atom q smarts
                parameter 1 [*:1] 0.0 # comment

            topology_procedure atom q sqm_am1bcc_elf10rdkit
                conformers 100
                filter_elf 10
                neutralize true
        
        model charge_transfer_bond
            topology_term c charge e float 0.0
            topology_procedure bond c smarts
                parameter 1 [#0:1]~[X2:2](~[*])~[*] -.1 # comment

        model charge_transfer_angle
            topology_term c charge e float 0.0
            topology_procedure angle c smarts
                parameter 1 [*:1]~[X2:2]([#0])~[*:3] -.1 # comment

        model charge_transfer_outofplane
            topology_term c charge e float 0.0
            topology_procedure outofplane c smarts
                parameter 1 [*:1]~[*:2](~[*:3])~[*:4] -.01 # comment

        model scale
            topology_term b scale12 none float 0.0
            topology_term a scale13 none float 0.0
            topology_term t scale14 none float 0.5
            topology_procedure torsion b a t smarts
                parameter 1 [*:1]~[*:2]~[*:3]~[*:4] 0.0 0.0 0.5 # comment

    chemical_model type=vdW name=""

        force type=lennard_jones_12_6_berkowitz_cutoff
            system_terms
                system_term id=c name=cutoff unit=A type=float default=10.0

            topology_terms
                topology_term id=e name=depth unit=kj/mol type=float
                topology_term id=s name=sigma unit=A type=float

            assignments:
                assignment type=smarts id=a topology=atom
                    <provides id=e/> 
                    <provides id=s/> 
                    <parameter id=1 smarts=[*:1] 
                        <value id=e>.001</value>
                        <value id=s>1.5</value>
                        <value id=e>.001</value>
                        <value id=s>1.5</value>
                    </parameter>
        force

        force scale
            topology_term b scale12 none float 0.0
            topology_term a scale13 none float 0.0
            topology_term t scale14 none float 0.5
            smarts a torsion b a t
                parameter 1 [*:1]~[*:2]~[*:3]~[*:4] 0.0 0.0 0.5 # comment


    chemical_model virtual_sites A

        transcoder t bond_to_angle
            hierarchy a
                parameter -1.5 180 [#6:1]~[#9:2]-[#0:3] ; #0 is a lone pair

        transcoder t angle_to_torsion
            hierarchy a
                parameter -1.3 120 180 [#1:1]-[#8:2]-[#1:3]-[#0:4] ; #0 is a lone pair

        transcoder t angle_to_tripyramid
            hierarchy a
                parameter -1.0 20 [#1:1]-[#8:2]([#0:4])([#8:5])-[#1:3] ; #0 is a lone pair

        force A angle_constrain
            hierarchy a
                parameter 180 [*:1]~[#9:2]-[#0:3] ; #0 is a lone pair

        force B bond_constrain
            hierarchy a
                parameter 1.5 [*]~[#9:1]-[#0:2] ; #0 is a lone pair

        force B bond_harmonic
            topology_term l length kJ/mol float
            topology_term k stiffness kJ/mol/A/A float

            hierarchy b bond k l
                parameter 1 1 1.5 [*:1]~[*:2] # comment

            hierarchy b bond k l
                parameter 1 1 1.5 [*:1]~[*:2] # comment

    hierarchy b angle angle_charge_shift
    parameter 1 1 0 [*:1]~[*:2]~[*:3] ; comment
"""
