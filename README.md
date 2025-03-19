# Implementing-802.11-Transmitter-and-Receiver

## Overview

This repository implements an OFDM (Orthogonal Frequency Division Multiplexing) communication system in MATLAB. It represents a simplified version of the WiFi (802.11) PHY layer, including packet construction at the transmitter and packet detection, synchronization, and decoding at the receiver.

## Features

Transmitter Implementation (OFDM_TX.m): Constructs and transmits OFDM packets, including preambles, pilots, and data symbols.

Receiver Implementation (OFDM_RX.m): Performs packet detection, synchronization, and decoding of received signals.

Modulation & Demodulation (mapping.m & demapper.m): Handles Forward Error Correction (FEC) and modulation/demodulation processes.

Channel Effects Simulation: Includes CFO (Carrier Frequency Offset), AWGN (Additive White Gaussian Noise), and other impairments.

## Files Description

OFDM_TX.m: Implements the transmitter, generating an OFDM waveform with preambles, IFFT processing, cyclic prefix addition, and interpolation.

OFDM_RX.m: Implements the receiver, handling packet detection, CFO estimation and correction, channel estimation, and payload processing.

mapping.m: Performs modulation mapping for data symbols.

demapper.m: Demodulates received signals and extracts data bits
