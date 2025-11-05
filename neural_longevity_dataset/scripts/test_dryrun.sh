#!/bin/bash

# Navigate to module root (parent of scripts/)
cd "$(dirname "$0")/.."

echo "🧬 Testing Longevity Pipeline (Dry-Run)"
echo ""

# Clean old data
echo "🧹 Cleaning old data..."
rm -rf longevity_dataset/
echo "✓ Clean"
echo ""

# Run dry-run
echo "🚀 Running dry-run..."
python3 neural_longevity_dataset.py --config configs/default.yaml --dry-run

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📊 Checking generated files:"
echo ""

if [ -f "longevity_dataset/samples_list.txt" ]; then
    echo "✓ samples_list.txt created"
    echo "  First 5 samples:"
    head -5 longevity_dataset/samples_list.txt | sed 's/^/    /'
    echo ""
else
    echo "✗ samples_list.txt NOT found"
fi

if [ -f "longevity_dataset/longevous_samples.txt" ]; then
    echo "✓ longevous_samples.txt created"
    echo "  First longevous sample:"
    head -1 longevity_dataset/longevous_samples.txt | sed 's/^/    /'
    echo ""
else
    echo "✗ longevous_samples.txt NOT found"
fi

if [ -f "longevity_dataset/central_points.json" ]; then
    echo "✓ central_points.json created"
    echo "  Number of points:"
    python3 -c "import json; data=json.load(open('longevity_dataset/central_points.json')); print(f'    {len(data)} central points')"
    echo ""
else
    echo "✗ central_points.json NOT found"
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "✨ Test complete!"
echo ""
echo "If everything is OK, you can:"
echo "  1. Configure API key in configs/default.yaml"
echo "  2. Run without --dry-run to process for real"
echo ""
echo "See: README.md for complete documentation"

