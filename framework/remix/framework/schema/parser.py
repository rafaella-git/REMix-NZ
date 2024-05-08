import json
import os
import re
from copy import deepcopy
from pathlib import Path

import yaml

from remix.framework import __gamscode__
from remix.framework import __versionhome__

# MANAGE USING A REMIX ENVIRONMENTAL VARIABLE?
MODEL_VERSION = "latest"
TARGET_DIRECTORY_INPUT_YAML = Path(__versionhome__).joinpath("schema").joinpath("input").joinpath("yaml")
TARGET_DIRECTORY_INPUT_JSON = Path(__versionhome__).joinpath("schema").joinpath("input").joinpath("json")
TARGET_DIRECTORY_OUTPUT_YAML = Path(__versionhome__).joinpath("schema").joinpath("output").joinpath("yaml")
TARGET_DIRECTORY_OUTPUT_JSON = Path(__versionhome__).joinpath("schema").joinpath("output").joinpath("json")
TARGETS = [TARGET_DIRECTORY_INPUT_YAML, TARGET_DIRECTORY_INPUT_JSON, TARGET_DIRECTORY_OUTPUT_YAML, TARGET_DIRECTORY_OUTPUT_JSON]
YEARTYPES = []  # ["years", "vintage", "yearsToCalc"]
KEYS = ["name", "title", "description", "constraints", "type", "units", "isAbout"]
REPEAT_SUFFIXES = ["", "_a"]
PROFILE_FIELDS = ["against", "along", "fixed", "lower", "upper"]
PROFILE_VALUE = {"name": "Value", "title": "Profile values", "type": "number"}
OEO_BASE_IRI = "http://openenergy-platform.org/ontology/oeo/"
MISSING_TERM_TRACKING = "https://github.com/OpenEnergyPlatform/ontology/issues"
BFO_BASE_IRI = "http://purl.obolibrary.org/obo/"


def main(as_resource=True):
    """Main schema parsing method.

    Parameters
    ----------
    as_resource : bool, optional
        If False, the output schemas will contain only OEP-compatible information, by default True.
    """
    # create version-specific folder structure
    os.makedirs(__versionhome__, exist_ok=True)
    file_list = [
        os.path.join(root, cur_file)
        for root, _, files in os.walk(__gamscode__)
        for cur_file in files
        if (".gms" in cur_file)
    ]
    [os.makedirs(t, exist_ok=True) for t in TARGETS]

    sets_file = next(f for f in file_list if "sets.gms" in f)
    with open(sets_file, "r") as f:
        lines = f.readlines()
    chunk = "".join(lines)

    set_data = get_all_sets(chunk)
    set_profile = get_profile_set(chunk)
    set_maps = get_all_maps(chunk)
    set_unions = get_all_unions(chunk)

    set_data.extend(set_unions)
    set_data.extend(set_profile)
    set_data = [expand_links(s) for s in set_data]
    set_schemas = [build_set_schema(s) for s in set_data]

    set_equalities = get_all_eqs(chunk)
    map_schemas = [build_map_schemas(m, set_data, set_equalities) for m in set_maps]
    # Add custom set for scenario linksModel and techs
    # set_data.append(["scenario", "Scenario", [{"name": "scenario", "path": "http://openenergy-platform.org/ontology/oeo/OEO_00000364"}],"set_scenarios", "Scenario indexes help to differentiate scenarios"])
    # set_data.append(["techs", "All technologies", [{"name": "energy transformation unit", "path": "http://openenergy-platform.org/ontology/oeo/OEO_00020102"}],"set_techs", "Concatenation of all technology sets"])
    # set_data.append(["linksModel", "Links Model", [{"name": "grid component link", "path": "http://openenergy-platform.org/ontology/oeo/OEO_00000255"}],"set_linksmodel", "Aggregated model links"])
    # set_data.append(["capType", "Capacity Type", [{"name": "grid component link", "path": "http://openenergy-platform.org/ontology/oeo/OEO_00000255"}],"set_linksmodel", "Aggregated model links"])
    for gams_file in file_list:
        # Extract data from the .gms files
        tables = extract_source_data(gams_file)
        outputs = extract_output_info(gams_file)
        if "sets" in gams_file:
            tables = [t for t in tables if "### map_" not in t]
        fields_fks = [parse_data(t) for t in tables]
        # NOT GENERAL == UGLY
        if "converter" in gams_file:
            ffk_coeffprofile = [field for field in fields_fks if field[0] == "converter_coefficientprofile"][0]
            idx = fields_fks.index(ffk_coeffprofile)
            ffk_coeffprofile= list(ffk_coeffprofile)
            ffk_coeffprofile[-1][-1] = ffk_coeffprofile[-1][-1].split(")")[0]
            ffk_coeffprofile[-2] = []
            fields_fks[idx] = tuple(ffk_coeffprofile)
        # END OF UGLYINESS
        schemas = [
            build_schemas(f, set_data, set_equalities, set_profile) for f in fields_fks
        ]

        write_schemas(schemas, as_resource,
                      target_json=TARGET_DIRECTORY_INPUT_JSON,
                      target_yaml=TARGET_DIRECTORY_INPUT_YAML)
        if len(outputs) > 0:
            output_fields = [parse_output_data(os) for os in outputs]
            out_schemas = [
            build_schemas(f, set_data, set_equalities, set_profile) for f in output_fields
        ]
            write_schemas(out_schemas, as_resource,
                      target_json=TARGET_DIRECTORY_OUTPUT_JSON,
                      target_yaml=TARGET_DIRECTORY_OUTPUT_YAML)

    write_schemas(map_schemas, as_resource,
                      target_json=TARGET_DIRECTORY_INPUT_JSON,
                      target_yaml=TARGET_DIRECTORY_INPUT_YAML)
    write_schemas(set_schemas, as_resource,
                      target_json=TARGET_DIRECTORY_INPUT_JSON,
                      target_yaml=TARGET_DIRECTORY_INPUT_YAML)


def expand_links(set_data):
    """Selects the ontology term from the set data tuple to pass it to the
    ontology expander.

    Parameters
    ----------
    set_data : tuple
        Tuple with set information

    Returns
    -------
    str
        Expanded ontology IRI
    """
    set_data_new = set_data
    set_data_new[2] = expand_ontology_link(set_data_new[2])
    return set_data_new


def expand_ontology_link(concept):
    """Expands ontology terms with their appropiate IRIs

    Parameters
    ----------
    concept : str
        Ontology identifier, only compatible with OEO and BFO by now.

    Returns
    -------
    str
        Ontology IRI of the term, missing terms get an issue tracker link.
    """
    path, name = concept.split(":")
    if "OEO_" in path:
        path = f"{OEO_BASE_IRI}{path}"
    if "BFO_" in path:
       path = f"{BFO_BASE_IRI}{path}"
    if "MISSING_TERM" in path:
        path = MISSING_TERM_TRACKING
        name = f"missing:{name}"
    return [{"name": name, "path":path}]


def write_schemas(schemas, as_resource=True,
                  target_json=TARGET_DIRECTORY_INPUT_JSON,
                  target_yaml=TARGET_DIRECTORY_INPUT_YAML):
    """Write schemas to .yaml and .json.

    Parameters
    ----------
    schemas : list
        List containing the schema name and respective data.
    as_resource : bool, optional
        If False, the output schemas will contain only OEP-compatible information, by default True.
    target_json : str
        Target directory of json schemas
    target_yaml : str
        Target directory of yaml schemas
    """
    for n, t, d, s in schemas:
        # Names should be always lower, according to tabular data resource spec.
        n = n.lower()
        resource = {
            "name": n,
            "title": t,
            "description": d,
            "schema": s
        }
        filepath = target_yaml / "{}.schema.yaml".format(n)
        output_file = resource if as_resource else s
        with open(filepath, 'w', encoding="utf-8") as output:
            yaml.dump(output_file, output, sort_keys=False, default_flow_style =False)
        filepath = target_json / "{}.schema.json".format(n)
        with open(filepath, 'w', encoding="utf-8") as output:
            json.dump(output_file, output, indent=4)


def extract_output_info(source_file: Path):
    """Open source_file and extract all the output file related data.

    Parameters
    ----------
    source_file : Path
        Path to the .gms source file to be read.

    Returns
    -------
    all_outputsx
        List of raw data of the output data.
    """
    with open(source_file, "r") as f:
        lines = f.readlines()
    chunk = "".join(lines)
    all_outputs = re.findall(r'(?s)(\*\s\/\/\sOUTPUT\s*:\s*.*?)(?=\"[ \t]*\;)', chunk, re.IGNORECASE)
    if len(all_outputs) == 0:
        all_outputs = re.findall(r'(?s)(\*\s\/\/\sOUTPUT\s*:\s*.*?)(?=\"\s*\n\d*)', chunk, re.IGNORECASE)
    return all_outputs

def extract_source_data(source_file: Path):
    """Open source_file and extract the markdown templates [markdown], and the raw parameter tables [tables].

    Parameters
    ----------
    source_file : Path
        Path to the .gms source file to be read.

    Returns
    -------
    all_tables
        List of raw data of the tables.
    """
    with open(source_file, "r") as f:
        lines = f.readlines()
    chunk = "".join(lines)
    # old_pattern:  (?s)((\*\s\/\/\s###\s\w+\n\*\s\/\/\sTitle.*?)(?=/;).*?)(?=\))\)
    all_tables = re.findall(r'(?s)((\*\*\s\/\/\sINPUT\s*:\s*.*?)(?=/;).*?)(?=\))\)', chunk, re.IGNORECASE)
    all_tables = [dom[0] for dom in all_tables]
    return all_tables


def get_fks(chunk):
    """Get foreign keys from the table data.

    Parameters
    ----------
    chunk : str
        Raw table data.

    Returns
    -------
    list
        foreign keys.
    """
    fks = chunk.split("\ntable")[1]
    fks = fks.split("(")[1].split(",")
    if ")" in fks[-1]:
        fks[-1] = fks[-1].split(")")[0]
    fks = [fk for fk in fks if "pc_" not in fk]
    return fks


def get_fields(chunk):
    """Get the fields from raw table data.

    Parameters
    ----------
    chunk : str
        Raw table data.

    Returns
    -------
    list
        List with field information of the table.
    """
    raw_list = chunk.split("\n")
    elements = [l.strip() for l in raw_list if len(l.strip()) != 0]
    elements = [e for e in elements if re.match(r"\w+\s+\"\w+", e)]
    elements = [re.sub(r"\s+\"", "|", e).strip(",").split("|") for e in elements]
    elements = [
        [sub.replace('"', "").strip(" ") for sub in e if len(sub) > 0] for e in elements
    ]

    fields = [dict(zip(KEYS, element)) for element in elements]
    for field in fields:
        if "isAbout" in field:
            field["isAbout"] = expand_ontology_link(field["isAbout"])
        if "type" not in field:
            field["type"] = "number"
        if ("constraints" in field) and (field["constraints"] != ""):
            constraints = field["constraints"].split(";")
            constraints = [c.split(":") for c in constraints]
            field["constraints"] = {
                c[0]: eval(c[1]) for c in constraints
            }
        else:
            field.pop("constraints", None)

        if ("constraints" in field) and (len(field["constraints"]) == 0):
            field.pop("constraints", None)

    return fields

def parse_data(chunk):
    """Transform the list of parameter strings into fields of the schema.

    Parameters
    ----------
    chunk : str
        Single table raw data

    Returns
    -------
    tuple
        name, title, description, fields and foreign_keys of this table.
    """
    # old_pattern: (?s)(\*\s\/\/\s\{table_.*?)(?=/;)
    table = re.findall(r"(?s)(\*\s\/\/\s###\s\w+.*\*\s*\/\/\s*Title.*?)(?=\/;)", chunk, re.IGNORECASE)[0]
    fields = get_fields(table)
    raw_list = table.split("\n")
    name = (
        raw_list[0]
        .replace("* //", "")
        .strip()
        .replace("{", "")
        .replace("}", "")
        .replace("### ", "")
    )
    title = (
        raw_list[1]
        .replace("* //", "")
        .strip()
        .replace("Title: ", "")
    )
    description = (
        raw_list[2]
        .replace("* //", "")
        .strip()
        .replace("Description: ", "")
    )
    foreign_keys = get_fks(chunk)
    return name, title, description, fields, foreign_keys

def parse_output_data(chunk):
    """Transform the list of parameter strings into fields of the schema.

    Parameters
    ----------
    chunk : str
        Single table raw data

    Returns
    -------
    tuple
        name, title, description, fields and foreign_keys of this table.
    """
    name = re.findall(r"(?s)(OUTPUT\:\s*(.*)?)((?=\s*\|))", chunk)[0][1].strip()
    title = re.findall(r"(?s)(\/\/\s*Title\:(.*)?)((?=\s*\n+\s*))", chunk)[0][1].strip().split("\n")[0]
    description = chunk.split("\"")[-1]
    isabout = expand_ontology_link(chunk.split("|")[-1].split("\n* //")[0].strip())
    fields = [{"name": "value", "title": title, "desription": description, "type": "number", "isAbout": isabout}]
    foreign_keys = re.findall(r"(?s)\(((.*)?)((?=\)))", chunk)[0][1].replace("%scenidx%", "scenario,").replace(" ", "").replace("\n", "").split(",")
    foreign_keys = [f.rstrip("_a").strip() for f in foreign_keys]
    return name, title, description, fields, foreign_keys


def build_schemas(name_fields_fks, set_data, equalities, profile_set):
    name, title, description, fields, fks = name_fields_fks
    has_profiletypes = (
        any([f["name"] in PROFILE_FIELDS for f in fields]) and "profile" in name.lower()
    )
    if has_profiletypes:
        fks.insert(-1, "profileTypes")
    fks_temp = [fk if fk not in equalities else equalities[fk] for fk in fks]
    set_maps = {s[0]: s for s in set_data if s[0] in fks_temp}
    for old, new in equalities.items():
        if old in fks:
            set_maps[old] = deepcopy(set_maps[new])
    # Left-right suffixes
    if any([len([k for k in fks if k == f]) > 1 for f in fks]):
        repeated = set([f for f in fks if len([k for k in fks if k == f]) > 1])
        for r in repeated:
            # indices = [i for i, x in enumerate(fks) if x == r]
            if len([k for k in fks if k == r]) > 2:
                raise ValueError
            else:
                new_fks = [r + REPEAT_SUFFIXES[i] for i in range(2)]
                for nf in new_fks:
                    set_maps[nf] = deepcopy(set_maps[r])
                fks = [f for f in fks if f != r]
                new_fks.extend(fks)
                fks = new_fks
    # Start-end suffixes
    if any("_start" in f for f in fks) and any("_end" in f for f in fks):
        start = [f for f in fks if "_start" in f][0]
        end = [f for f in fks if "_end" in f][0]
        if start.rstrip("_start") != end.rstrip("_end"):
            raise ValueError
        else:
            extended_map = {s[0]: s for s in set_data if s[0] in [start.rstrip("_start")]}
            set_maps[start] = deepcopy(extended_map[start.rstrip("_start")])
            set_maps[end] = deepcopy(extended_map[start.rstrip("_start") ])
    fields = (
        [f for f in fields if f["name"] not in PROFILE_FIELDS]
        if has_profiletypes
        else fields
    )
    if "profile" in name.lower():
        fields.extend([PROFILE_VALUE])

    foreign_key_fields = []
    schema_fks = []

    for fk in fks:
        foreign_key_fields += [{
            "name": fk,
            "title": set_maps[fk][1],
            "type": "integer" if fk in YEARTYPES else "string",
            "isAbout": set_maps[fk][2]
        }]
        schema_fks += [{
            "fields": [fk],
            "reference": {
                "fields": [set_maps[fk][0]],
                "resource": set_maps[fk][3]
            },
        }]

    foreign_key_fields.extend(fields)
    schema = {"fields": foreign_key_fields, "foreignKeys": schema_fks}
    return name, title, description, schema


def get_all_sets(chunk):
    """Extract all sets from the sets.gms source file.

    Parameters
    ----------
    chunk : str
        Contents of the sets.gms file.

    Returns
    -------
    list
        Raw set data.
    """
    # old pattern: (?s)(\*\*\s+\/\/\s+SET.*?)(?=.csv)
    sets = re.findall(r"(?s)(\*\*\s+\/\/\s+SET.*?)(?=\/)", chunk, re.IGNORECASE)
    sets_data = [s.split(r"SET:")[-1].strip().split("|") for s in sets]
    for set_d in sets_data:
        set_d[-1], set_description = set_d[-1].split(".csv")
        set_description = set_description.split("\"")[1]
        set_d.append(set_description)
    sets_data = [[sub.strip() for sub in s] for s in sets_data]
    return sets_data


def get_profile_set(chunk):
    """Extract all profile sets from the sets.gms source file.

    Parameters
    ----------
    chunk : str
        Contents of the sets.gms file.

    Returns
    -------
    list
        Raw set data.
    """
    sets = re.findall(r"(?s)(\*\*\s+\/\/\s+PROFILE.*?)(?=.csv)", chunk, re.IGNORECASE)
    profile_set = [[z.strip() for z in s.split(r"PROFILE:")[-1].strip().split("|")] for s in sets]
    profile_set
    profile_set[0].append("Labels for profile levels")
    return profile_set


def build_set_schema(set_data):
    name, title, concept, extension, description = set_data
    type = "string" if name not in YEARTYPES else "number"
    schema = {
        "fields": [
            {
                "name": name,
                "title": title,
                "type": type,
                "constraints": {"unique": True},
                "isAbout": concept
            }
        ],
        "primaryKey": [name],
    }
    return extension, title, description, schema


def get_all_maps(chunk):
    """Extract all maps from the sets.gms source file.

    Parameters
    ----------
    chunk : str
        Contents of the sets.gms file.

    Returns
    -------
    list
        Raw set data.
    """
    maps = re.findall(r"(?s)(\*\*\s+\/\/\s+MAP.*?)(?=\nset)", chunk, re.IGNORECASE)
    map_data = [s.split(r"MAP:")[-1].strip().split("|") for s in maps]
    for mapd_d in map_data:
        mapd_d[-1], map_info = mapd_d[-1].split(".csv")
        map_info = map_info.split("\n")
        map_description = [i for i in map_info if "Description: " in i][0].split("Description: ")[-1]
        map_title = [i for i in map_info if "Title: " in i][0].split("Title: ")[-1]
        mapd_d.append(map_description)
        mapd_d.append(map_title)
    map_data = [[sub.strip() for sub in m] for m in map_data]
    return map_data


def get_all_eqs(chunk):
    """Extract all equalities from the sets.gms source file.

    Parameters
    ----------
    chunk : str
        Contents of the sets.gms file.

    Returns
    -------
    list
        Raw set data.
    """
    sames = re.findall(r"(?s)(\*\*\s+\/\/\s+SAME.*?)(?=\n)", chunk, re.IGNORECASE)
    sames_data = [s.split(r"SAME:")[-1].strip().split("|") for s in sames]
    sames_data = [[sub.strip() for sub in s] for s in sames_data]
    sames_data = {s[0]: s[1] for s in sames_data}
    return sames_data


def get_all_unions(chunk):
    """Extract all unions from the sets.gms source file.

    Parameters
    ----------
    chunk : str
        Contents of the sets.gms file.

    Returns
    -------
    list
        Raw set data.
    """
    # old pattern: (?s)(\*\*\s+\/\/\s+UNION.*?)(?=.csv)
    unions = re.findall(r"(?s)(\*\*\s+\/\/\s+UNION.*?)(?=\/)", chunk, re.IGNORECASE)
    union_data = [s.split(r"UNION:")[-1].strip().split("|") for s in unions]
    for union_d in union_data:
        union_d[-1], union_description = union_d[-1].split(".csv")
        union_description = union_description.split("\"")[1]
        union_d.append(union_description)
    union_data = [[sub.strip() for sub in u] for u in union_data]
    return union_data


def build_map_schemas(map_, sets, equalities):
    title = map_[-1]
    description = map_[-2]
    map_ = map_[:-2]
    map_ = [mp if mp not in equalities else equalities[mp] for mp in map_]
    set_maps = {m[0]: m[1:] for m in sets if m[0] in map_[0:-1]}
    fields = [
        {"name": m, "title": set_maps[m][0], "type": "string"} for m in map_[0:-1]
    ]
    foreignKeys = [
        {"fields": [m], "reference": {"fields": [m], "resource": set_maps[m][2]}}
        for m in map_[0:-1]
    ]
    name = map_[-1]
    schema = {"fields": fields, "foreignKeys": foreignKeys}
    return name , title, description, schema


if __name__ == "__main__":
    main()
