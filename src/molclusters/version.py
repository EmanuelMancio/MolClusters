# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Set version of package."""

import importlib.metadata

version = importlib.metadata.version("molclusters")
__version__ = version
