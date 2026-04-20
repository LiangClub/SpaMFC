"""
Feature Fusion Module

Supports three fusion methods:
- adaptive: Use Random Forest to calculate feature importance
- fixed: Use pre-defined fixed weights
- concat: Direct concatenation without weighting

Users can choose fusion method via configuration:
--fusion-method adaptive/fixed/concat
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from typing import Dict, Optional, List
import warnings


VALID_FUSION_METHODS = ["adaptive", "fixed", "concat"]


class FeatureFusion:
    """Feature fusion class supporting multiple methods"""
    
    def __init__(self, fusion_method: str = "adaptive"):
        """
        Initialize feature fusion
        
        Parameters:
            fusion_method: "adaptive", "fixed", or "concat"
        """
        if fusion_method not in VALID_FUSION_METHODS:
            warnings.warn(f"Unknown fusion method '{fusion_method}', using 'adaptive'")
            fusion_method = "adaptive"
        
        self.fusion_method = fusion_method
    
    def fuse(
        self,
        feature_dfs: Dict[str, pd.DataFrame],
        fixed_weights: Optional[Dict[str, float]] = None
    ) -> pd.DataFrame:
        """
        Fuse features using specified method
        
        Parameters:
            feature_dfs: Dictionary of feature DataFrames
                        {"spatial": df, "cnv": df, "expression": df, "niche": df}
            fixed_weights: Dictionary of fixed weights (used when fusion_method="fixed")
        
        Returns:
            Fused feature DataFrame
        """
        if self.fusion_method == "concat":
            return self._concat_features(feature_dfs)
        elif self.fusion_method == "fixed":
            weights = self._calculate_fixed_weights(feature_dfs, fixed_weights)
            return self._weighted_fusion(feature_dfs, weights)
        else:
            weights = self._calculate_adaptive_weights(feature_dfs)
            return self._weighted_fusion(feature_dfs, weights)
    
    def _concat_features(
        self,
        feature_dfs: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """Direct concatenation without weighting"""
        active_features = [k for k, v in feature_dfs.items() if v is not None]
        
        if len(active_features) == 0:
            return None
        
        fused_parts = []
        
        for feature_name in active_features:
            df = feature_dfs[feature_name]
            
            scaler = StandardScaler()
            scaled = scaler.fit_transform(df.values)
            
            scaled_df = pd.DataFrame(
                scaled,
                index=df.index,
                columns=[f"{feature_name}_{col}" for col in df.columns]
            )
            
            fused_parts.append(scaled_df)
        
        fused_df = pd.concat(fused_parts, axis=1)
        fused_df = fused_df.fillna(0.0)
        
        return fused_df
    
    def _weighted_fusion(
        self,
        feature_dfs: Dict[str, pd.DataFrame],
        weights: Dict[str, float]
    ) -> pd.DataFrame:
        """Weighted feature fusion"""
        active_features = [k for k, v in feature_dfs.items() if v is not None]
        
        if len(active_features) == 0:
            return None
        
        fused_parts = []
        
        for feature_name in active_features:
            df = feature_dfs[feature_name]
            
            scaler = MinMaxScaler()
            scaled = scaler.fit_transform(df.values)
            
            scaled_df = pd.DataFrame(
                scaled,
                index=df.index,
                columns=[f"{feature_name}_{col}" for col in df.columns]
            )
            
            weight = weights.get(feature_name, 1.0 / len(active_features))
            weighted = scaled_df * weight
            
            fused_parts.append(weighted)
        
        fused_df = pd.concat(fused_parts, axis=1)
        fused_df = fused_df.fillna(0.0)
        
        return fused_df
    
    def _calculate_fixed_weights(
        self,
        feature_dfs: Dict[str, pd.DataFrame],
        fixed_weights: Optional[Dict[str, float]] = None
    ) -> Dict[str, float]:
        """Calculate fixed weights"""
        active_features = [k for k, v in feature_dfs.items() if v is not None]
        
        if fixed_weights is None:
            default_weights = {
                "spatial": 0.3,
                "cnv": 0.3,
                "expression": 0.4,
                "niche": 0.2
            }
            fixed_weights = default_weights
        
        weights = {}
        total_weight = 0.0
        
        for feature_name in active_features:
            if feature_name in fixed_weights:
                weights[feature_name] = fixed_weights[feature_name]
                total_weight += fixed_weights[feature_name]
        
        if total_weight > 0:
            weights = {k: v / total_weight for k, v in weights.items()}
        else:
            weights = {k: 1.0 / len(active_features) for k in active_features}
        
        return weights
    
    def _calculate_adaptive_weights(
        self,
        feature_dfs: Dict[str, pd.DataFrame]
    ) -> Dict[str, float]:
        """Calculate adaptive weights using Random Forest"""
        active_features = [k for k, v in feature_dfs.items() if v is not None]
        
        if len(active_features) == 0:
            return {}
        
        if len(active_features) == 1:
            return {active_features[0]: 1.0}
        
        fused_df = pd.concat(
            [feature_dfs[k] for k in active_features],
            axis=1
        )
        
        fused_df = fused_df.fillna(0.0)
        
        pseudo_labels = self._generate_pseudo_labels(feature_dfs, fused_df)
        
        try:
            rf = RandomForestClassifier(
                n_estimators=200,
                random_state=0,
                max_depth=10
            )
            rf.fit(fused_df.values, pseudo_labels)
            
            feature_importances = rf.feature_importances_
            
            weights = {}
            start_idx = 0
            
            for feature_name in active_features:
                feature_dim = feature_dfs[feature_name].shape[1]
                end_idx = start_idx + feature_dim
                
                importance = np.mean(feature_importances[start_idx:end_idx])
                weights[feature_name] = importance
                
                start_idx = end_idx
            
            total_importance = sum(weights.values())
            if total_importance > 0:
                weights = {k: v / total_importance for k, v in weights.items()}
            else:
                weights = {k: 1.0 / len(active_features) for k in active_features}
            
            return weights
            
        except Exception as e:
            warnings.warn(f"Random Forest failed: {e}, using equal weights")
            return {k: 1.0 / len(active_features) for k in active_features}
    
    def _generate_pseudo_labels(
        self,
        feature_dfs: Dict[str, pd.DataFrame],
        fused_df: pd.DataFrame
    ) -> np.ndarray:
        """Generate pseudo labels for Random Forest training"""
        n_samples = len(fused_df)
        
        if "spatial" in feature_dfs and feature_dfs["spatial"] is not None:
            spatial_df = feature_dfs["spatial"]
            
            if "distance_to_center" in spatial_df.columns:
                distance = spatial_df["distance_to_center"]
                if distance.std() > 0:
                    pseudo_labels = (distance > distance.median()).astype(int).values
                    return pseudo_labels
            
            if spatial_df.shape[1] > 0:
                first_col = spatial_df.iloc[:, 0]
                if first_col.std() > 0:
                    pseudo_labels = (first_col > first_col.median()).astype(int).values
                    return pseudo_labels
        
        pseudo_labels = np.random.randint(0, 2, n_samples)
        return pseudo_labels
    
    def get_weights(
        self,
        feature_dfs: Dict[str, pd.DataFrame],
        fixed_weights: Optional[Dict[str, float]] = None
    ) -> Dict[str, float]:
        """Get weights for current fusion method"""
        if self.fusion_method == "concat":
            return {k: 1.0 for k in feature_dfs.keys() if feature_dfs[k] is not None}
        elif self.fusion_method == "fixed":
            return self._calculate_fixed_weights(feature_dfs, fixed_weights)
        else:
            return self._calculate_adaptive_weights(feature_dfs)
    
    def get_weight_summary(self, weights: Dict[str, float]) -> str:
        """Get weight summary string"""
        summary = f"Fusion method: {self.fusion_method}\n"
        summary += "Feature weights: "
        for feature_name, weight in weights.items():
            summary += f"{feature_name}={weight:.3f}, "
        return summary.rstrip(", ")


class AdaptiveWeightCalculator:
    """Legacy class for backward compatibility"""
    
    def __init__(self, weight_method: str = "adaptive"):
        self.fusion = FeatureFusion(fusion_method=weight_method)
    
    def calculate(
        self,
        feature_dfs: Dict[str, pd.DataFrame],
        fixed_weights: Optional[Dict[str, float]] = None
    ) -> Dict[str, float]:
        return self.fusion.get_weights(feature_dfs, fixed_weights)
    
    def normalize_features(
        self,
        feature_dfs: Dict[str, pd.DataFrame],
        weights: Dict[str, float]
    ) -> pd.DataFrame:
        return self.fusion._weighted_fusion(feature_dfs, weights)
    
    def get_weight_summary(self, weights: Dict[str, float]) -> str:
        return self.fusion.get_weight_summary(weights)