# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Set version of package."""

import importlib.metadata

version = importlib.metadata.version("mdrhconstant")
__version__ = version
