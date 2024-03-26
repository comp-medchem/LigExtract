{
  "query": {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_nonpolymer_entity_instance_container_identifiers.comp_id",
          "operator": "exact_match",
          "value": "LIGAND"
        }
      },
      {
        "logical_operator": "or",
        "type": "group",
        "nodes": [
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "attribute": "rcsb_nonpolymer_instance_feature_summary.comp_id",
              "operator": "exists",
              "negation": true
            }
          },
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "attribute": "rcsb_nonpolymer_instance_feature_summary.type",
              "operator": "exact_match",
              "value": "HAS_COVALENT_LINKAGE",
              "negation": true
            }
          },
          {
            "logical_operator": "and",
            "type": "group",
            "nodes": [
              {
                "type": "terminal",
                "service": "text",
                "parameters": {
                  "attribute": "rcsb_nonpolymer_instance_feature_summary.comp_id",
                  "operator": "exact_match",
                  "value": "LIGAND",
                  "negation": true
                }
              },
              {
                "type": "terminal",
                "service": "text",
                "parameters": {
                  "attribute": "rcsb_nonpolymer_instance_feature_summary.type",
                  "operator": "exact_match",
                  "value": "HAS_COVALENT_LINKAGE",
                  "negation": true
                }
              }
            ]
          }
        ]
      }
    ]
  },
  "return_type": "non_polymer_entity",
  "request_options": {
    "results_verbosity": "compact"
  }
}
