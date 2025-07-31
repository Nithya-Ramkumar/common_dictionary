#!/usr/bin/env python3
"""
Test script for ChemSpider API with polymer names
"""

import requests
import json
import time
from typing import Dict, Any, List

class ChemSpiderTester:
    def __init__(self, api_key: str):
        self.api_key = api_key
        self.base_url = "https://www.chemspider.com"
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })
    
    def search_by_name(self, name: str) -> List[int]:
        """Search for compounds by name"""
        url = f"{self.base_url}/Search.asmx/SimpleSearch"
        params = {
            'query': name,
            'token': self.api_key
        }
        
        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse XML response
            import xml.etree.ElementTree as ET
            root = ET.fromstring(response.text)
            
            # Extract CSIDs
            csids = []
            for csid_elem in root.findall('.//int'):
                csids.append(int(csid_elem.text))
            
            print(f"Search for '{name}': Found {len(csids)} CSIDs")
            return csids
            
        except Exception as e:
            print(f"Error searching for '{name}': {e}")
            return []
    
    def get_compound_info(self, csid: int) -> Dict[str, Any]:
        """Get basic compound information"""
        url = f"{self.base_url}/Compound.asmx/GetCompoundInfo"
        params = {
            'CSID': csid,
            'token': self.api_key
        }
        
        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse XML response
            import xml.etree.ElementTree as ET
            root = ET.fromstring(response.text)
            
            info = {}
            for elem in root:
                tag = elem.tag.replace('{http://www.chemspider.com}', '')
                info[tag] = elem.text
            
            return info
            
        except Exception as e:
            print(f"Error getting compound info for CSID {csid}: {e}")
            return {}
    
    def get_synonyms(self, csid: int) -> List[str]:
        """Get synonyms for a compound"""
        url = f"{self.base_url}/Compound.asmx/GetSynonyms"
        params = {
            'CSID': csid,
            'token': self.api_key
        }
        
        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse XML response
            import xml.etree.ElementTree as ET
            root = ET.fromstring(response.text)
            
            synonyms = []
            for synonym_elem in root.findall('.//string'):
                synonyms.append(synonym_elem.text)
            
            return synonyms
            
        except Exception as e:
            print(f"Error getting synonyms for CSID {csid}: {e}")
            return []
    
    def get_properties(self, csid: int) -> Dict[str, Any]:
        """Get molecular properties"""
        url = f"{self.base_url}/Compound.asmx/GetProperties"
        params = {
            'CSID': csid,
            'token': self.api_key
        }
        
        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse XML response
            import xml.etree.ElementTree as ET
            root = ET.fromstring(response.text)
            
            properties = {}
            for prop_elem in root.findall('.//Property'):
                name = prop_elem.find('Name')
                value = prop_elem.find('Value')
                if name is not None and value is not None:
                    properties[name.text] = value.text
            
            return properties
            
        except Exception as e:
            print(f"Error getting properties for CSID {csid}: {e}")
            return {}

def test_polymer_names():
    """Test ChemSpider with polymer names"""
    
    # You'll need to provide your API key
    api_key = input("Enter your ChemSpider API key: ").strip()
    
    if not api_key:
        print("No API key provided. Exiting.")
        return
    
    tester = ChemSpiderTester(api_key)
    
    # Test polymer names
    polymer_names = ["PVDF", "PBI", "PI", "Polyvinylidene fluoride", "Polybenzimidazole", "Polyimide"]
    
    for name in polymer_names:
        print(f"\n{'='*50}")
        print(f"Testing: {name}")
        print(f"{'='*50}")
        
        # Search by name
        csids = tester.search_by_name(name)
        
        if csids:
            # Get info for first few results
            for i, csid in enumerate(csids[:3]):  # Limit to first 3 results
                print(f"\n--- Result {i+1} (CSID: {csid}) ---")
                
                # Get compound info
                info = tester.get_compound_info(csid)
                if info:
                    print("Compound Info:")
                    for key, value in info.items():
                        if value:
                            print(f"  {key}: {value}")
                
                # Get synonyms
                synonyms = tester.get_synonyms(csid)
                if synonyms:
                    print(f"Synonyms ({len(synonyms)}):")
                    for synonym in synonyms[:5]:  # Show first 5 synonyms
                        print(f"  - {synonym}")
                
                # Get properties
                properties = tester.get_properties(csid)
                if properties:
                    print("Properties:")
                    for key, value in properties.items():
                        if value:
                            print(f"  {key}: {value}")
                
                # Add delay to avoid rate limiting
                time.sleep(1)
        else:
            print("No results found")
        
        # Add delay between searches
        time.sleep(2)

if __name__ == "__main__":
    test_polymer_names() 