"""picbh tests."""

import picbh


def test_stub() -> None:
    """Stub test to ensure the test suite runs."""
    print(picbh.__version__)  # noqa: T201


def test__greet() -> None:
    """Test the greet function."""
    assert picbh.greet("World") == "Hello, World!"


def test__greet_jim() -> None:
    """Test the greet_jim function."""
    assert picbh.greet_jim() == "Hello, Jim!"
