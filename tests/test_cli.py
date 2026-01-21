import pytest
import sys
import os
from unittest.mock import patch, MagicMock
from pdb_cli import cli

def test_format_value():
    # Test simple dict formatting
    val = {"A": 1, "B": 2}
    formatted = cli.format_value(val)
    assert "A: 1" in formatted
    assert "B: 2" in formatted

    # Test list formatting
    val_list = ["Item 1", "Item 2"]
    formatted_list = cli.format_value(val_list)
    assert "Item 1" in formatted_list
    assert "Item 2" in formatted_list

@patch("pdb_cli.pdb_summary.summarize")
def test_main_summary(mock_summarize):
    # Mock return value
    mock_summarize.return_value = {"Test Feature": "Test Value"}
    
    with patch.object(sys, "argv", ["pdb-cli", "test.pdb", "--chains"]):
        cli.main()
        
    mock_summarize.assert_called_once()
    args, kwargs = mock_summarize.call_args
    assert args[0] == "test.pdb"
    assert kwargs["chains"] is True

@patch("requests.get")
def test_download_pdb(mock_get, tmp_path):
    # Mock successful download
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.content = b"fake pdb content"
    mock_get.return_value = mock_response

    outdir = tmp_path
    cli.download_pdb("1abc", outdir=str(outdir))
    
    expected_file = outdir / "1abc.pdb"
    assert expected_file.exists()
    assert expected_file.read_bytes() == b"fake pdb content"

@patch("requests.get")
def test_download_alphafold(mock_get, tmp_path):
    # Mock successful download
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.content = b"fake af content"
    mock_get.return_value = mock_response

    outdir = tmp_path
    cli.download_alphafold("P12345", outdir=str(outdir))
    
    expected_file = outdir / "AF-P12345-F1-model_v6.pdb"
    assert expected_file.exists()
    assert expected_file.read_bytes() == b"fake af content"
