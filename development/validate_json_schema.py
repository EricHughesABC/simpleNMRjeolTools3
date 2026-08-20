#!/usr/bin/env python3
"""
NMR JSON Data Validator

This program validates JSON files containing NMR spectroscopy data against
a predefined JSON schema. It provides detailed error reporting and can be
used both as a standalone script and as an imported module.

Requirements:
    pip install jsonschema

Usage:
    python nmr_validator.py input_file.json
    python nmr_validator.py input_file.json --schema custom_schema.json
"""

import json
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import jsonschema
from jsonschema import validate, ValidationError, Draft7Validator


# The NMR data schema definition
NMR_SCHEMA = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "type": "object",
    "title": "NMR Spectroscopy Data Schema",
    "description": "Schema for NMR spectroscopy experimental data and assignments",
    "additionalProperties": {
        "anyOf": [
            {
                "type": "object",
                "properties": {
                    "datatype": {"type": "string", "const": "nmrspectrum"},
                    "origin": {"type": "string"},
                    "type": {"type": "string", "enum": ["1D", "2D"]},
                    "subtype": {"type": "string"},
                    "experimenttype": {"type": "string"},
                    "experiment": {"type": "string"},
                    "class": {"type": "string"},
                    "spectype": {"type": "string"},
                    "pulsesequence": {"type": "string"},
                    "intrument": {"type": "string"},
                    "probe": {"type": "string"},
                    "datafilename": {"type": "string"},
                    "nucleus": {
                        "oneOf": [
                            {"type": "string"},
                            {"type": "array", "items": {"type": "string"}}
                        ]
                    },
                    "specfrequency": {
                        "oneOf": [
                            {"type": "number"},
                            {"type": "array", "items": {"type": "number"}},
                            {"type": "string"}
                        ]
                    },
                    "temperature": { 
                        "oneOf": [
                            {"type": "string"},
                            {"type": "number"}
                        ]
                    },
                    "peaks": {
                        "type": "object",
                        "properties": {
                            "datatype": {"type": "string", "const": "peaks"},
                            "count": {"type": "integer"},
                            "data": {
                                "type": "object",
                                "patternProperties": {
                                    "^[0-9]+$": {
                                        "type": "object",
                                        "properties": {
                                            "intensity": {"type": "number"},
                                            "type": {"type": "integer"},
                                            "annotation": {"type": "string"},
                                            "delta1": {"type": "number"},
                                            "delta2": {"type": "number"}
                                        },
                                        "required": ["intensity", "type", "annotation", "delta1", "delta2"]
                                    }
                                },
                                "additionalProperties": False
                            }
                        },
                        "required": ["datatype", "count", "data"]
                    },
                    "integrals": {
                        "type": "object",
                        "properties": {
                            "datatype": {"type": "string", "const": "integrals"},
                            "count": {"type": "integer"},
                            "normValue": {"type": "number"},
                            "data": {"type": "object"}
                        },
                        "required": ["datatype", "count", "normValue", "data"]
                    },
                    "multiplets": {
                        "type": "object",
                        "properties": {
                            "datatype": {"type": "string", "const": "multiplets"},
                            "count": {"type": "integer"},
                            "normValue": {"type": "number"},
                            "data": {"type": "object"}
                        },
                        "required": ["datatype", "count", "normValue", "data"]
                    },
                    "filename": {"type": "string"}
                },
                "required": [
                    "datatype", "peaks", "integrals", "multiplets"
                ],
                "additionalProperties": True
            },
            {
                "type": "object",
                "properties": {
                    "datatype": {"type": "string"},
                    "count": {"type": "integer"},
                    "data": {"type": "object"}
                },
                "required": ["datatype", "count", "data"],
                "additionalProperties": False
            }
        ]
    }
}


class NMRJSONValidator:
    """Validator class for NMR JSON data files."""
    
    def __init__(self, schema: Optional[Dict] = None):
        """
        Initialize the validator with a schema.
        
        Args:
            schema: JSON schema dictionary. If None, uses the default NMR schema.
        """
        self.schema = schema or NMR_SCHEMA
        self.validator = Draft7Validator(self.schema)
    
    def load_json_file(self, filepath: str) -> Tuple[Optional[Dict], Optional[str]]:
        """
        Load JSON data from a file.
        
        Args:
            filepath: Path to the JSON file
            
        Returns:
            Tuple of (data, error_message). If successful, error_message is None.
        """
        try:
            with open(filepath, 'r', encoding='utf-8') as file:
                data = json.load(file)
            return data, None
        except FileNotFoundError:
            return None, f"File not found: {filepath}"
        except json.JSONDecodeError as e:
            return None, f"Invalid JSON format: {e}"
        except Exception as e:
            return None, f"Error reading file: {e}"
    
    def validate_data(self, data: Dict) -> Tuple[bool, List[str]]:
        """
        Validate JSON data against the NMR schema.
        
        Args:
            data: JSON data dictionary
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        
        try:
            validate(instance=data, schema=self.schema)
            return True, []
        except ValidationError as e:
            errors.append(self._format_validation_error(e))
            
            # Collect all validation errors, not just the first one
            for error in self.validator.iter_errors(data):
                formatted_error = self._format_validation_error(error)
                if formatted_error not in errors:
                    errors.append(formatted_error)
            
            return False, errors
    
    def _format_validation_error(self, error: ValidationError) -> str:
        """
        Format a validation error into a readable string.
        
        Args:
            error: ValidationError instance
            
        Returns:
            Formatted error message string
        """
        path = " -> ".join(str(p) for p in error.absolute_path) if error.absolute_path else "root"
        return f"Error at '{path}': {error.message}"
    
    def validate_file(self, filepath: str) -> Tuple[bool, List[str], Optional[Dict]]:
        """
        Validate a JSON file against the NMR schema.
        
        Args:
            filepath: Path to the JSON file
            
        Returns:
            Tuple of (is_valid, error_messages, data)
        """
        data, load_error = self.load_json_file(filepath)
        
        if load_error:
            return False, [load_error], None
        
        is_valid, validation_errors = self.validate_data(data)
        return is_valid, validation_errors, data
    
    def get_data_summary(self, data: Dict) -> Dict[str, any]:
        """
        Generate a summary of the NMR data.
        
        Args:
            data: Validated JSON data dictionary
            
        Returns:
            Dictionary containing data summary
        """
        summary = {
            "total_properties": len(data),
            "spectra_count": 0,
            "total_atoms": 0,
            "carbon_atoms": 0,
            "experiment_types": []
        }
        
        # Count spectra by looking for objects with datatype = nmrspectrum
        for key, value in data.items():
            if isinstance(value, dict) and value.get("datatype") == "nmrspectrum":
                summary["spectra_count"] += 1
        
        # Get atom counts
        if "allAtomsInfo" in data:
            summary["total_atoms"] = data["allAtomsInfo"].get("count", 0)
        
        if "carbonAtomsInfo" in data:
            summary["carbon_atoms"] = data["carbonAtomsInfo"].get("count", 0)
        
        # Get experiment types
        if "exptIdentifiers" in data:
            exp_data = data["exptIdentifiers"].get("data", {})
            summary["experiment_types"] = list(exp_data.values())
        
        return summary


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Validate NMR JSON data files against the NMR schema",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python nmr_validator.py data.json
    python nmr_validator.py data.json --schema custom_schema.json
    python nmr_validator.py data.json --verbose
        """
    )
    
    parser.add_argument(
        "input_file",
        help="Path to the JSON file to validate"
    )
    
    parser.add_argument(
        "--schema",
        help="Path to custom JSON schema file (optional)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed data summary after validation"
    )
    
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Only show validation result (pass/fail)"
    )
    
    args = parser.parse_args()
    
    # Load custom schema if provided
    schema = None
    if args.schema:
        try:
            with open(args.schema, 'r', encoding='utf-8') as f:
                schema = json.load(f)
            if not args.quiet:
                print(f"Using custom schema: {args.schema}")
        except Exception as e:
            print(f"Error loading custom schema: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Initialize validator
    validator = NMRJSONValidator(schema)
    
    # Validate the file
    is_valid, errors, data = validator.validate_file(args.input_file)
    
    if args.quiet:
        # Quiet mode: just exit with appropriate code
        sys.exit(0 if is_valid else 1)
    
    # Print results
    print(f"Validating: {args.input_file}")
    print("=" * 50)
    
    if is_valid:
        print("VALIDATION PASSED")
        print("The JSON file is valid according to the NMR schema.")
        
        if args.verbose and data:
            print("\nData Summary:")
            print("-" * 20)
            summary = validator.get_data_summary(data)
            for key, value in summary.items():
                print(f"{key.replace('_', ' ').title()}: {value}")
    else:
        print("VALIDATION FAILED")
        print(f"Found {len(errors)} error(s):")
        print()
        
        for i, error in enumerate(errors, 1):
            print(f"{i}. {error}")
    
    sys.exit(0 if is_valid else 1)


if __name__ == "__main__":
    main()