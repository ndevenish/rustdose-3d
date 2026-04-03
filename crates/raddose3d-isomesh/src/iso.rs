/// Iso-level selection mirroring `computeIsoLevels(maxDose)` from `web/index.html:684-698`.
pub struct IsoLevel {
    pub value: f32,
    pub color: [f32; 3],
    pub opacity: f32,
}

// Named colour constants matching the JS
const LIGHTBLUE: [f32; 3] = [0.678, 0.847, 0.902];
const ROYALBLUE: [f32; 3] = [0.255, 0.412, 0.882];
const RED: [f32; 3] = [1.0, 0.0, 0.0];

pub fn compute_iso_levels(max_dose: f32) -> Vec<IsoLevel> {
    if max_dose <= 0.0 {
        return vec![];
    }
    if max_dose > 30.0 {
        vec![
            IsoLevel {
                value: 0.1,
                color: LIGHTBLUE,
                opacity: 0.2,
            },
            IsoLevel {
                value: 20.0,
                color: ROYALBLUE,
                opacity: 0.5,
            },
            IsoLevel {
                value: 30.0,
                color: RED,
                opacity: 0.9,
            },
        ]
    } else {
        vec![
            IsoLevel {
                value: (max_dose * 0.10_f32).max(1e-6),
                color: LIGHTBLUE,
                opacity: 0.2,
            },
            IsoLevel {
                value: max_dose * 0.50,
                color: ROYALBLUE,
                opacity: 0.5,
            },
            IsoLevel {
                value: max_dose * 0.85,
                color: RED,
                opacity: 0.9,
            },
        ]
    }
}
