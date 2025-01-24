{
  "query": {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_nonpolymer_instance_annotation.comp_id",
          "operator": "exact_match",
          "value": "LIGAND"
        }
      },
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_nonpolymer_instance_annotation.type",
          "operator": "exact_match",
          "value": "HAS_NO_COVALENT_LINKAGE"
        }
      }
    ]
  },
  "return_type": "non_polymer_entity",
  "request_options": {
    "results_verbosity": "compact"
  }
}
