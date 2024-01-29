use std::{
    f32::consts::PI,
    ops::Range,
};

use plotters::prelude::*;
use rand::random;
use sogi_pll::{PllConfig, PllResult, SogiPll};

const SAMPLING_TIME: f32 = 1.0 / 12000.0;
const END: f32 = PI / 24.0;
const RANGE: Range<f32> = 0.0..END;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let omega_n = 2.0 * PI * 50.0;

    let config = PllConfig {
        sample_time: SAMPLING_TIME,
        sogi_k: 1.0,
        pi_proportional_gain: 178.0,
        pi_integral_gain: 0.0001,
        omega_zero: omega_n,
    };

    let mut pll = SogiPll::new(config);

    let x = RANGE.step(SAMPLING_TIME);

    let y_cos: Vec<f32> = x
        .values()
        .map(|x| {
            (2.0 * PI * 50.0 * x + random::<f32>() * 1.5 - 1.23).sin()
                + (random::<f32>() - 0.5) * 0.30
        })
        .collect();

    let pll_out: Vec<PllResult> = y_cos.clone().into_iter().map(|x| pll.update(x)).collect();

    let theta: Vec<f32> = pll_out.iter().map(|x| x.theta).collect();
    let rms: Vec<f32> = pll_out.iter().map(|x| x.v_rms()).collect();

    let reconstructed_wave = theta.clone().into_iter().map(|y| (y).cos());

    let root = BitMapBackend::new("SOGI PLL Demo.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0.0..END, -2f32..2f32)?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(x.values().zip(y_cos), RED))?
        .label("Input")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));

    chart
        .draw_series(LineSeries::new(x.values().zip(theta), GREEN))?
        .label("Theta")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GREEN));

    chart
        .draw_series(LineSeries::new(x.values().zip(reconstructed_wave), BLACK))?
        .label("Reconstructed")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK));

    chart
        .draw_series(LineSeries::new(x.values().zip(rms), BLUE))?
        .label("RMS")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}
