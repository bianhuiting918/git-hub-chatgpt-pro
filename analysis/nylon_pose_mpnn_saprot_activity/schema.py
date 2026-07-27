ENDPOINTS = [
    "pa66_l1_color_span_px",
    "pa66_l2_color_span_px",
    "pa6_max_color_span_px",
    "l1_fraction",
]
PRIMARY_PREDICTORS = ["pose_strict_nac", "mpnn_vanilla_mean", "saprot_full_mean"]
SENSITIVITY_PREDICTORS = ["pose_loose_gate", "pose_weight5_loose", "pose_weight10_loose", "mpnn_soluble_mean"]
ENDPOINT_LABELS = {
    "pa66_l1_color_span_px": "PA66 L1 raster color-span proxy (px)",
    "pa66_l2_color_span_px": "PA66 L2 raster color-span proxy (px)",
    "pa6_max_color_span_px": "PA6 maximum raster color-span proxy (px)",
    "l1_fraction": "PA66 L1 fraction of L1+L2 raster spans",
}
