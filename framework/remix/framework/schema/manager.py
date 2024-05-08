import json
import logging
import os
import re
from glob import glob
from pathlib import Path
from pathlib import PurePath
from typing import Union

import jsonschema
import requests

from remix.framework import __remixhome__
from remix.framework import __version__
from remix.framework import __versionhome__
from remix.framework.schema.parser import main as parse

LATEST = Path(__versionhome__).joinpath("schema").joinpath("input").joinpath("json")
METADATA_EXISTS = LATEST.exists()
METADATA_ISEMPTY = not any(LATEST.iterdir()) if METADATA_EXISTS else True
LOCAL_PROFILE = Path(__remixhome__).joinpath("specifications").joinpath("frictionless").joinpath("table-schema.json")
PROFILE = "https://specs.frictionlessdata.io/schemas/table-schema.json"


class Manager:
    """REMix data schema manager. This class has methods to get the schemas from their input files."""

    def __init__(self):
        # create .remix user folder structure
        if not os.path.exists(__remixhome__):
            os.makedirs(__remixhome__, exist_ok=True)
        # Check schema consistency by init
        self.check_latest()

    def check_latest(self):
        """This method checks for consistency of the generated schema files."""
        if (not METADATA_EXISTS) or METADATA_ISEMPTY:
            parse()
        # validate
        internal_check = LATEST.exists()
        if internal_check:
            if not LOCAL_PROFILE.parent.is_dir():
                LOCAL_PROFILE.parent.mkdir(exist_ok=True, parents=True)
            if not LOCAL_PROFILE.exists():
                response = requests.get(PROFILE)
                if response.status_code == 200:
                    validation_profile = response.json()
                    with open(LOCAL_PROFILE, "w") as fp:
                        json.dump(validation_profile, fp)
                else:
                    logging.INFO(
                        "The metadata profile could not be downloaded, skipping validation"
                    )
                    return
            else:
                with open(LOCAL_PROFILE, "r") as fp:
                    validation_profile = json.load(fp)
            validator_class = jsonschema.validators.validator_for(validation_profile)
            validator = validator_class(validation_profile)
            try:
                schemas = self.get_schemas("latest")
            except KeyError:
                logging.INFO(
                    "The local schemas are invalid, deleting and downloading again"
                )
                parse()
            errors = {}
            for schema in schemas:
                try:
                    validator.validate(schema)
                except jsonschema.ValidationError:
                    errors[schema] = []
                    for e in validator.iter_errors(schema):
                        errors[schema].append(e)
            if len(errors.keys()) != 0:
                wrong_schemas = list(errors.keys())
                raise jsonschema.ValidationError(
                    f"There was errors in: {wrong_schemas}"
                )

    def get_resources(self, version="latest"):
        """Get the model schemas as they are written by the parser.

        This includes keys that are not valid in the context of a tabular data schema but are necessary for a resource
        schema.

        Parameters
        ----------
        version : str, optional
            Version of the model resource schemas, by default "latest".

        Returns
        -------
        dict
            Dictionary with the resource schemas of the model.

        Raises
        ------
        KeyError
            If an invalid version is given.
        """
        available_versions = list(Path(__remixhome__).iterdir())
        available_versions = [s for s in available_versions if (bool(re.search(r"\d*\.\d*\.\d*", str(s)))) or (f"{__version__}" in str(s))]
        schema_dict = {PurePath(s).parts[-1]: s for s in available_versions}
        latest_version = LATEST.as_posix()

        if version == "latest":
            schema_paths = self.find_schema_paths(latest_version)
        else:
            if version in schema_dict:
                schema_paths = self.find_schema_paths(schema_dict[version])
            else:
                raise KeyError(f"Schemas version {version} is not valid")

        schemas = {
            PurePath(s).parts[-1].split(".")[0]: self.load_schema(s)
            for s in schema_paths
        }
        schemas["version"] = (
            version if version != "latest" else PurePath(latest_version).parts[-1]
        )
        return schemas

    def get_schemas(self, version: str = "latest") -> dict:
        """Returns schemas valid in the context of a tabular data schema.

        Parameters
        ----------
        version : str, optional
            Version of the model schemas, by default "latest".

        Returns
        -------
        dict
            Dictionary with the schemas of the model.
        """
        resources = self.get_resources(version=version)
        schemas = {k: v["schema"] for k, v in resources.items() if k != "version"}
        return schemas

    def get_tabular_package_template(self, version: str = "latest") -> dict:
        """Get the given schema version as a template for building tabular data packages.

        Parameters
        ----------
        version : str, optional
            version of the model schemas, by default "latest".

        Returns
        -------
        dict
            Dictionary with the tabular data packages.
        """
        resources = self.get_resources(version=version)
        resource_array = [v for k, v in resources.items() if k != "version"]
        ver = resources["version"]
        package = {
            "name": f"REMix data v{ver}",
            "profile": "tabular-data-package",
            "resources": resource_array,
        }
        return package

    @staticmethod
    def find_schema_paths(path: str) -> dict:
        """Helper method to find schemas in a given directory.

        Parameters
        ----------
        path : str
            Path of the directory containing the schemas.

        Returns
        -------
        list
            List of the schemas contained in the given directory.
        """
        schema_paths = glob(os.path.join(path, "*.schema.json"))
        return schema_paths

    @staticmethod
    def load_schema(path) -> dict:
        """Helper method to load schemas from its path.

        Parameters
        ----------
        path : str
            Path to the schema.

        Returns
        -------
        dict
            Dictionary with a single schema information.
        """
        with open(os.path.join(path), "rb") as f:
            schema = json.loads(f.read())
        return schema

    @staticmethod
    def calculate_value(string: str) -> Union[str, int]:
        """Calculate value of the schema semantic version.

        Parameters
        ----------
        string : str
            Version name.

        Returns
        -------
        str, int
            Version number as integer or "latest" if version is "latest".
        """
        vals = [int(st) for st in string.split(".") if st != "latest"]
        return (
            vals[0] * 100 + vals[1] * 10 + vals[2] * 1 if string != "latest" else string
        )
