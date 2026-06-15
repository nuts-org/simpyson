from __future__ import annotations

import re
from abc import ABC, abstractmethod
from typing import Any


class PulseSequenceTemplate(ABC):
    """Base class for pulse sequence templates"""

    def __init__(self, **kwargs):
        """
        Initialize pulse sequence template.

        Args:
            **kwargs: Parameters to override defaults (without variable_ prefix)
        """

        self.parameters = self.get_default_parameters()

        prefixed_params = {}
        for key, value in kwargs.items():
            prefixed_key = f"variable_{key}" if not key.startswith('variable_') else key
            prefixed_params[prefixed_key] = value

        self.parameters.update(prefixed_params)

        # Validate parameters
        self.validate_parameters()

    @abstractmethod
    def get_default_parameters(self) -> dict[str, Any]:
        """Return default parameters for this pulse sequence (with variable_ prefix)"""

    @abstractmethod
    def get_required_parameters(self) -> set[str]:
        """Return set of required parameter names (with variable_ prefix)"""

    def validate_parameters(self):
        """Validate that all required parameters are present"""
        missing = self.get_required_parameters() - set(self.parameters.keys())
        if missing:
            missing_clean = {param.replace('variable_', '') for param in missing}
            raise ValueError(f"Missing required parameters: {', '.join(missing_clean)}")

    def update_parameters(self, **kwargs):
        """
        Update parameters and re-validate.

        Args:
            **kwargs: Parameters to update (without variable_ prefix)
        """

        prefixed_params = {}
        for key, value in kwargs.items():
            prefixed_key = f"variable_{key}" if not key.startswith('variable_') else key
            prefixed_params[prefixed_key] = value

        self.parameters.update(prefixed_params)
        self.validate_parameters()

class NoPulse(PulseSequenceTemplate):
    """
    No pulse sequence - direct acquisition.

    Parameters:
        tsw (float): Dwell time in microseconds. Default: '1e6/sw'
    """

    def get_default_parameters(self) -> dict[str, Any]:
        return {
            'variable_tsw': '1e6/sw',
            'variable_offset': 0.0,
            'variable_num_channels': 1,
        }

    def get_required_parameters(self) -> set[str]:
        return {'variable_tsw'}

    @property
    def description(self) -> str:
        return "No pulse, direct acquisition"

    def generate_code(self) -> str:
        offset_val = self.parameters.get('variable_offset', 0.0)
        n_chan = int(self.parameters.get('variable_num_channels', 1))
        if offset_val == 0:
            offset_line = ""
        elif n_chan <= 1:
            offset_line = "    offset $par(offset)\n"
        else:
            extra = " 0" * (n_chan - 1)
            offset_line = f"    offset $par(offset){extra}\n"
        return f"""
proc pulseq {{}} {{
    global par
{offset_line}    acq_block {{
        delay $par(tsw)
    }}
}}
"""

class Pulse90(PulseSequenceTemplate):
    """
    Single 90° pulse on 1H.

    On multi-channel spin systems, SIMPSON's ``pulse`` command requires an
    rf/phase pair per channel; extra channels are padded with ``0 0``.

    Parameters:
        pH (float): Pulse length in microseconds. Default: 5.0
        plH (float): Pulse power in Hz. Default: 50000
        phH (str): Pulse phase. Default: '90'
        tsw (float): Dwell time in microseconds. Default: '1e6/sw'
    """

    def get_default_parameters(self) -> dict[str, Any]:
        return {
            'variable_pH': 5.0,
            'variable_plH': 50000,
            'variable_phH': '90',
            'variable_tsw': '1e6/sw',
            'variable_num_channels': 1,
        }

    def get_required_parameters(self) -> set[str]:
        return {'variable_pH', 'variable_plH', 'variable_phH', 'variable_tsw'}

    @property
    def description(self) -> str:
        return "Single 90° pulse on 1H"

    def generate_code(self) -> str:
        # SIMPSON requires one "rf phase" pair per channel on the pulse line
        n_chan = int(self.parameters.get('variable_num_channels', 1))
        extra = " 0 0" * (n_chan - 1)
        return f"""
proc pulseq {{}} {{
    global par
    pulse $par(pH) $par(plH) $par(phH){extra}
    acq_block {{
        delay $par(tsw)
    }}
}}
"""

class CPMAS(PulseSequenceTemplate):
    """
    Cross-polarization magic angle spinning sequence.

    The initial 1H magnetization is set via ``start_operator`` (e.g. ``I1x``)
    in the SIMPSON par block, so no explicit 90° pulse is needed here.

    Note: ``sw`` does not need to be commensurate with
    ``spin_rate * gamma_angles``; SIMPSON samples ``acq_block`` at ``1/sw``
    independently of the block duration (verified empirically).

    Parameters:
        pcp (float): Contact pulse length in μs. Default: 1000
        plHcp (float): 1H contact pulse power in Hz. Default: 70000
        phHcp (str): 1H contact pulse phase. Default: '0'
        plCcp (float): 13C contact pulse power in Hz. Default: 69000
        phCcp (str): 13C contact pulse phase. Default: '0'
        dw (str): Dwell time expression. Default: '1e6/spin_rate/gamma_angles'
    """

    def __init__(self, **kwargs):
        # Populated by SimpCalc from the actual spinsys before generate_code() is called.
        self.turnoff_interactions: list[str] = []
        super().__init__(**kwargs)

    def get_default_parameters(self) -> dict[str, Any]:
        return {
            'variable_pcp': 1000,
            'variable_plHcp': 70000,
            'variable_phHcp': '0',
            'variable_plCcp': 69000,
            'variable_phCcp': '0',
            'variable_dw': '1e6/spin_rate/gamma_angles'
        }

    def get_required_parameters(self) -> set[str]:
        return {
            'variable_pcp', 'variable_plHcp', 'variable_phHcp',
            'variable_plCcp', 'variable_phCcp', 'variable_dw'
        }

    @property
    def description(self) -> str:
        return "Cross-polarization magic angle spinning"

    def generate_code(self) -> str:
        turnoff_line = (
            f"    turnoff {' '.join(self.turnoff_interactions)}\n"
            if self.turnoff_interactions else ""
        )
        return (
            "\nproc pulseq {} {\n"
            "    global par\n"
            "    reset\n"
            "    pulse $par(pcp) $par(plHcp) $par(phHcp) $par(plCcp) $par(phCcp)\n"
            + turnoff_line
            + "    acq_block {\n"
            "        delay $par(dw)\n"
            "    }\n"
            "}\n"
        )

pulseq_templates = {
    'no_pulse': NoPulse,
    'pulse_90': Pulse90,
    'cp_mas': CPMAS,
}

def get_template(name: str, **kwargs) -> PulseSequenceTemplate:
    """
    Get a pulse sequence template by name.

    Args:
        name: Template name ('no_pulse', 'pulse_90', 'cp_mas')
        **kwargs: Parameters to override defaults

    Returns:
        PulseSequenceTemplate: Configured template instance
    """
    if name not in pulseq_templates:
        available = ', '.join(pulseq_templates.keys())
        raise ValueError(f"Unknown template '{name}'. Available: {available}")

    return pulseq_templates[name](**kwargs)

class CustomPulseSequence(PulseSequenceTemplate):
    """Wrapper for custom string pulse sequences"""

    def __init__(self, code: str, **kwargs):
        self.code = code
        # Filter kwargs to only include parameters that appear in the code

        required = self.get_required_parameters()
        filtered_kwargs = {}
        for k, v in kwargs.items():
            # Accept both already-prefixed keys and unprefixed keys
            if k in required or f"variable_{k}" in required:
                filtered_kwargs[k] = v

        super().__init__(**filtered_kwargs)

    def get_default_parameters(self) -> dict[str, Any]:
        return {}

    def get_required_parameters(self) -> set[str]:
        # Extract parameter names from the code using regex
        params = set()
        pattern = r'\$par\(([^)]+)\)'
        matches = re.findall(pattern, self.code)
        for match in matches:
            params.add(f"variable_{match}")
        return params

    @property
    def description(self) -> str:
        return "Custom pulse sequence"

    def generate_code(self) -> str:
        if self.code.strip().startswith("proc pulseq"):
            return self.code

        return f"""
proc pulseq {{}} {{
    global par
{self.code}
}}
"""
