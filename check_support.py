"""Checks that remain active under Python optimization."""

def require(condition, message="verification failed"):
    if not condition:
        raise AssertionError(message)
