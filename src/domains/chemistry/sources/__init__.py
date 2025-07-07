from .pubchem_source import PubChemSource
from .rdkit_source import RDKitSource
from .source_factory import SourceFactory
 
SourceFactory.register_source_type("pubchem", PubChemSource)
SourceFactory.register_source_type("rdkit", RDKitSource) 