from pydantic.v1 import BaseModel
from biochatter.api_agent.base.formatters import format_as_rest_call, format_as_python_call


class TestParams(BaseModel):
    """Test parameter model for REST calls."""

    base_url: str
    endpoint: str
    param1: str | None = None
    param2: int | None = None
    question_uuid: str | None = None


class TestMethodParams(BaseModel):
    """Test parameter model for Python api_name calls."""

    api_name: str
    param1: str | None = None
    param2: int | None = None
    question_uuid: str | None = None


def test_format_as_rest_call_basic():
    params = TestParams(
        base_url="http://api.example.com",
        endpoint="/v1/test",
        param1="value1",
        param2=42,
    )
    expected = "http://api.example.com/v1/test?param1=value1&param2=42"
    assert format_as_rest_call(params) == expected


def test_format_as_rest_call_with_trailing_slash():
    params = TestParams(
        base_url="http://api.example.com/",
        endpoint="/v1/test/",
        param1="value1",
    )
    expected = "http://api.example.com/v1/test?param1=value1"
    assert format_as_rest_call(params) == expected


def test_format_as_rest_call_no_params():
    params = TestParams(
        base_url="http://api.example.com",
        endpoint="/v1/test",
    )
    expected = "http://api.example.com/v1/test?"
    assert format_as_rest_call(params) == expected


def test_format_as_python_call_with_module():
    params = TestMethodParams(
        api_name="sqrt",
        param1="16",
        param2=42,
    )
    expected = "sqrt(param1='16', param2=42)"
    assert format_as_python_call(params) == expected


def test_format_as_python_call_without_module():
    params = TestMethodParams(
        api_name="calculate",
        param1="test",
        param2=123,
    )
    expected = "calculate(param1='test', param2=123)"
    assert format_as_python_call(params) == expected


def test_format_as_python_call_no_params():
    params = TestMethodParams(
        api_name="empty_function",
    )
    expected = "empty_function()"
    assert format_as_python_call(params) == expected


def test_format_as_python_call_ignores_question_uuid():
    params = TestMethodParams(
        api_name="test",
        param1="value",
        question_uuid="123-456",
    )
    expected = "test(param1='value')"
    assert format_as_python_call(params) == expected
