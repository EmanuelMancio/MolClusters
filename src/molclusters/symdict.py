# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the `SymmetricDict` class."""

from collections.abc import Hashable, Iterable, MutableMapping
from typing import TypeVar, overload

_KT = TypeVar("_KT", bound=Hashable)
_VT = TypeVar("_VT")


class SymmetricDict(MutableMapping[_KT, _VT]):
    """
    Symmetrical dictionary structure.

    A dictionary-like data structure storing values for key pairs symmetrically,
    with single storage per unique pair.

    Parameters
    ----------
    None

    Notes
    -----
    - Access via tuples: x[a, b] or x[b, a] gives the same value.
    - Single-key access: x[a] returns a temporary dictionary of connected keys.
    - MutableMapping methods (.get, .keys, .items, etc.) are fully supported.
    """

    _data: dict[tuple[_KT, _KT], _VT]

    def __init__(self) -> None:
        """Initialize an empty SingleStoreDict."""
        self._data = {}

    def _normalize(self, a: _KT, b: _KT) -> tuple[_KT, _KT]:
        """
        Normalize a key pair for symmetric storage.

        Parameters
        ----------
        a : _KT
            First key.
        b : _KT
            Second key.

        Returns
        -------
        tuple[_KT, _KT]
            Ordered tuple (min(a,b), max(a,b)) for unique storage.
        """
        return (a, b) if a <= b else (b, a)

    @overload
    def __getitem__(self, key: _KT) -> dict[_KT, _VT]: ...

    @overload
    def __getitem__(self, key: tuple[_KT, _KT]) -> _VT: ...

    def __getitem__(self, key: _KT | tuple[_KT, _KT]) -> dict[_KT, _VT] | _VT:
        """
        Retrieve a value or a temporary dictionary of connected keys.

        Parameters
        ----------
        key : _KT or tuple[_KT, _KT]
            Single key for temporary dict view or tuple for direct value access.

        Returns
        -------
        _VT or dict[_KT, _VT]
            Value for a tuple key, or temporary dict of neighbors for a single key.
        """
        if isinstance(key, tuple):
            k1, k2 = key
            return self._data[self._normalize(k1, k2)]
        else:
            return {
                other: v
                for (k1, k2), v in self._data.items()
                if k1 == key or k2 == key
                for other in [(k2 if k1 == key else k1)]
            }

    @overload
    def __setitem__(self, key: _KT, value: dict[_KT, _VT]) -> None: ...
    @overload
    def __setitem__(self, key: tuple[_KT, _KT], value: _VT) -> None: ...

    def __setitem__(
        self, key: _KT | tuple[_KT, _KT], value: dict[_KT, _VT] | _VT
    ) -> None:
        """
        Assign a value to a key pair.

        Parameters
        ----------
        key : _KT or tuple[_KT, _KT]
            Single key (with dict value) or tuple key.
        value : _VT or dict[_KT, _VT]
            Value to store or dict of values for multiple connections.

        Raises
        ------
        ValueError
            If single key assignment is not a dictionary.
        """
        if isinstance(key, tuple):
            self._data[self._normalize(*key)] = value  # single pair
        elif isinstance(value, dict):
            for k, v in value.items():
                self._data[self._normalize(key, k)] = v
        else:
            raise ValueError("Single-key assignment must be a dict")

    def __delitem__(self, key: _KT | tuple[_KT, _KT]) -> None:
        """
        Delete a key or all pairs containing a key.

        Parameters
        ----------
        key : _KT or tuple[_KT, _KT]
            Single key (deletes all related pairs) or tuple key.
        """
        if isinstance(key, tuple):
            del self._data[self._normalize(*key)]
        else:
            for k in [pair for pair in self._data if key in pair]:
                del self._data[k]

    def __iter__(self) -> Iterable[tuple[_KT, _KT]]:
        """Iterate over all unique key pairs (normalized order only).

        Returns
        -------
        Iterable
        """
        return iter(self._data)

    def __len__(self) -> int:
        """Return the number of unique keys.

        Returns
        -------
        int
        """
        return len(self._data)

    def __repr__(self) -> str:
        """Return a string representation of the internal storage.

        Returns
        -------
        str
        """
        return repr(self._data)

    def all_keys(self) -> set[_KT]:
        """
        Return all individual keys stored in the dictionary.

        Returns
        -------
        set[_KT]
            Unique keys from all pairs.
        """
        return {k for pair in self._data for k in pair}
